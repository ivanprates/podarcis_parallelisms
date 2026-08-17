## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, March 2026.
## The goals of this R script are:
## To perform analyses of phenotypic trait variation, correlation, and trajectories in Podarcis lizards.
## To compose manuscript figures.

## 1. Getting ready ----

## Packages:
library(geomorph)
library(paletteer)
library(patchwork)
library(psych)
library(reshape2)
library(MCMCglmm)
library(ape)
library(coda)
library(tidyverse)
library(svglite)

## Path to working directory: IP's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(paste0(path, "morphological_analyses/"))

## Create folders to store outputs:
dir.create("data/")
dir.create("plots/")

## Coloration and morphometric data for all measured specimens (including some not sequenced):
idf <- read.csv(file = "data/Podarcis_greenbrowns_LizWiz_outputs_sample_info_2025-06-25.csv", header = TRUE, na.strings = "NA")

## If needed, keeping only genetic samples (30 per species) + P. muralis (also 30, for balanced sampling).
## Read list of sequenced samples + P. muralis:
gs <- read.csv(file = paste0(path, "sample_information/genetic_samples/samples_DNA_extractions_all_in_2024-02_plus_selected_muralis.csv"), header = TRUE)

## Keep only morphological data for sequenced specimens + P. muralis:
idf <- idf[idf$unique_identifier %in% unique(gs$unique_identifier), ]

## Exclude potential juveniles (if present)?
## Not needed if subsetting only to sequenced individuals.
juv <- unique(idf$unique_identifier[grepl("subadult|juvenile", idf$Notes) & !is.na(idf$Notes)]) ; juv
idf <- idf[idf$unique_identifier %in% setdiff(idf$unique_identifier, juv), ]

## Another criterion for juvenile exclusion -- small size:
## Adjust for each species if needed (e.g., P. liolepis is a fairly small species).
#ms <- 45
#ms <- 50
#small <- unique(idf$unique_identifier[idf$svl < ms & !is.na(idf$svl)]) ; small
#idf <- idf[idf$unique_identifier %in% setdiff(idf$unique_identifier, small), ]

## If needed, remove females present in the sample.
## IP focused on males only, but JA and RGR have female data also.
## Not needed if subsetting only to sequenced individuals.
idf <- idf[idf$sex == "M" & !is.na(idf$sex), ]

## IP's focal populations:
f <- c("muralis_AZ", "muralis_PL",
       "filfolensis_CI", "filfolensis_CO",
       "pityusensis_PUJ", "pityusensis_TRU",
       "liolepis_MUR", "liolepis_GOD",
       "vaucheri_OU", "vaucheri_TZ",
       "cretensis_THE", "cretensis_EIS",
       "gaigeae_LIN", "gaigeae_AT")

## Keep only focal populations:
idf <- idf[idf$species_abbpop %in% f, ]

## Set factor levels to arrange downstream plot elements by species and green-brown pair:
idf$species_abbpop <- factor(x = idf$species_abbpop, levels = f)

## Melt data (i.e., reformat to long format):
dm <- melt(idf, id.vars = c("unique_identifier", "species", "species_abbpop", "group", 
                            "svl", "headlength", "headlength_by_SVL", "headlength_resid", "slice_label"))

## List layers (i.e., combinations of body part and coloration slices):
sort(unique(dm$slice_label))

## Traits to keep from these layers:
traits <- c("ML_ratio", "u", "s", "m", "l", "lum", "invlum", "area_tot", "area_rel", "h.theta")

## Keep them:
dm <- dm[dm$variable %in% traits, ]

## Add labels for the traits from each slice:
dm$trait <- paste(dm$slice_label, dm$variable, sep = "_")

## Traits to keep for downstream analyses:
traits_selected <- c("dorsum_total_color_h.theta",
                     "dorsum_total_melanin_invlum",
                     "right_OVS_UV_u",
                     "right_OVS_UV_area_tot",
                     "venter_total_melanin_area_rel",
                     #"venter_total_melanin_area_tot",
                     "svl",
                     #"headlength")
                     "headlength_by_SVL")

## Keep only these traits:
ds <- dm[dm$trait %in% traits_selected, ]
ds <- ds[setdiff(names(ds), c("slice_label", "variable"))]

## Back to wide format:
dp <- pivot_wider(data = ds,
                  names_from = trait,
                  values_from = value,
                  id_cols = c(unique_identifier, species, species_abbpop, group, 
                              svl, headlength_by_SVL, headlength_resid, headlength))

## Square root of areas:
#dp$right_OVS_UV_area_tot <- sqrt(as.numeric(dp$right_OVS_UV_area_tot))
#dp$venter_total_melanin_area_tot <- sqrt(as.numeric(dp$venter_total_melanin_area_tot))

## Trait values as numeric:
dp[traits_selected] <- as.data.frame(sapply(dp[traits_selected], as.numeric))

## How many individuals (including those with incomplete data) from each species in dataset?
table(dp$species)

## Remove individuals with NA for any trait:
dp <- na.omit(dp)

## After removing individuals with NA for any trait, how many individuals from each species in dataset?
table(dp$species)

## Save:
#write.csv(dp, file = "data/Podarcis_greenbrowns_selected_traits_2026-03-15.csv", row.names = FALSE)

## 3. MCMCglmm ----

## Labels:
taxa <- rev(c("muralis", "filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"))
taxon_labels <- rev(c("P. muralis", "P. filfolensis", "P. pityusensis", "P. vaucheri", "P. liolepis", "P. cretensis", "P. gaigeae"))
trait_labels <- c("Dorsal hue", "Dorsal darkness", "UV-spots reflectance", 
                  "UV-spots area", "Ventral melanized area", 
                  "Snout-vent length", "Relative head length")

## Tree:
tree <- read.tree(paste0(path, "gen_vs_morph_analyses/tree/Yang_et_al_2021_wall_lizards_time-tree.nex"))

## Ladderize:
tree <- ladderize(tree)

## Extract focal taxa:
ptaxa <- c("cretensis", "filfolensis", "gaigeae", "liolepis", "muralis_IT", "pityusensis", "vaucheri_SS")
tree <- keep.tip(tree, ptaxa)
tree$tip.label <- gsub(x = tree$tip.label, pattern = "_.+", replacement = "") ## Renaming tips.
tree <- phytools::force.ultrametric(tree, method = "extend")
#plot(tree)

## Build inverse:
Ainv <- inverseA(tree)$Ainv

## Trait data:
dat <- dp[c("species", "group", traits_selected)]
dat <- as.data.frame(dat)

## Scale traits:
dat <- dat %>% mutate(across(all_of(traits_selected), ~ as.numeric(scale(.))))

## Factors:
dat$species <- factor(dat$species, levels = tree$tip.label)
dat$group <- factor(dat$group)

## MCMCglmm Model ----
# model_p <- MCMCglmm(
#   
#   ## Traits:
#   cbind(dorsum_total_color_h.theta,
#         dorsum_total_melanin_invlum,
#         right_OVS_UV_u, 
#         right_OVS_UV_area_tot,
#         venter_total_melanin_area_rel,
#         svl,
#         headlength_by_SVL) ~ 
#     
#   ## Intercepts: Vary by trait.
#   trait - 1 +
#   
#   ## Effects: Group and species.
#   trait:species +
#   trait:group:species,
#   
#   ## Random effect:
#   ## Phylogenetic covariance among species means.
#   random = ~ us(trait):species,
#   ginverse = list(species = Ainv),
#   
#   ## Residual covariance:
#   rcov   = ~ us(trait):units, ## Covariance among traits.
#   
#   ## Priors:
#   prior = list(
#     G = list(G1 = list(V = diag(7), nu = 8)),
#     R = list(V = diag(7), nu = 8)),
#   
#   ## Other params:
#   family = rep("gaussian", 7),
#   data = dat,
#   pr = TRUE,
#   nitt = 500000,
#   burnin = 100000,
#   thin = 100,
#   verbose = TRUE)
# 
## Check convergence:
#plot(model_1$Sol)  ## Fixed effects traces.
#plot(model_1$VCV)  ## Variance component traces.
#autocorr(model_1$Sol)
#autocorr(model_1$VCV)
# 
# ## Phylogenetic structure:
# View(summary(model_p)$Gcovariances)
# View(apply(model_p$VCV, 2, mean))
# 
# ## Extract effects:
# effects <- summary(model_p)$solutions
# 
# ## Edit:
# effects <- as.data.frame(effects)
# effects <- round(x = effects, digits = 3)
# effects$model_labels <- row.names(effects)
# row.names(effects) <- NULL
# 
# ## Extract within species group effects:
# GB_effects <- effects[grep(x = effects$model_labels, pattern = "groupGreen"), ]
# GB_effects <- GB_effects %>% separate(model_labels, into = c("trait", "species", "group"), sep = ":")
# GB_effects$trait <- gsub(x = GB_effects$trait, pattern = "trait", replacement = "")
# GB_effects$species <- gsub(x = GB_effects$species, pattern = "species", replacement = "")
# #View(GB_effects)
# 
# ## Reorder, relabel by defining factors:
# GB_effects$species <- factor(x = GB_effects$species, levels = taxa, labels = taxon_labels)
# GB_effects$trait <- factor(x = GB_effects$trait, levels = traits_selected, labels = trait_labels)
# 
# ## Save:
# write.csv(x = GB_effects, file = "data/MCMCglm_group_effects.csv", row.names = FALSE)

## If reading from a previous analysis:
GB_effects <- read.csv(file = "data/MCMCglm_group_effects.csv", header = TRUE)

## Factors:
GB_effects$species <- factor(x = GB_effects$species, levels = taxon_labels)
GB_effects$trait <- factor(x = GB_effects$trait, levels = trait_labels)

## Add categories:
GB_effects <- GB_effects %>%
  mutate(sig = case_when(
    l.95..CI > 0 ~ "Up in green",
    u.95..CI < 0 ~ "Up in brown",
    TRUE ~ "Not significant"))

## Add a dummy row so that muralis appears alone in first row (put legend there later):
GB_effects$species <- factor(GB_effects$species, levels = c("P. muralis", "", "P. filfolensis", "P. pityusensis","P. vaucheri","P. liolepis","P. cretensis","P. gaigeae"))

## Factor levels:
GB_effects$trait <- factor(GB_effects$trait, levels = rev(trait_labels))

## Plot forest:
f2a <- ggplot(GB_effects, aes(x = post.mean, y = trait, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = l.95..CI, xmax = u.95..CI),
                height = 0.2, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c(
    "Up in green" = "#31B565",
    "Up in brown" = "#9B5300",
    "Not significant" = "grey40")) +
  facet_wrap(~species, ncol = 2, drop = FALSE) +
  labs(x = "Posterior effect size", y = NULL, colour = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    strip.text = element_text(face = "italic", size = 12, colour = "black"),
    strip.background = element_rect(fill = "gray90", colour = "transparent"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.775, 0.915),
    panel.spacing = unit(0.75, "lines"))
  
## Check:
f2a

## 4. Testing trait correlations per species ----

## Function: correlation analyses:
run_co <- function(species) { 

  ## Testing:
  #species <- "cretensis"
  
  ## Status:
  print(paste0("Now processing ", species))
  
  ## Data for focal species:
  vd <- dp[dp$species == species, ]
  
  ## Subset columns:
  vd <- vd[setdiff(names(vd), c("unique_identifier", "species", "species_abbpop", "group"))]
  
  ## Subset traits:
  vd <- vd[traits_selected]
  
  ## As numeric:
  vd <- as.data.frame(sapply(vd, as.numeric))

  ## Test for trait correlations:
  res <- corr.test(vd, use = "pairwise", method = "pearson")
  
  ## Extract correlation coefficients and p-values:
  tc_c <- res$r
  tc_p <- res$p

  ## Long format:
  tc_c <- melt(tc_c, varnames = c("trait_1", "trait_2"), value.name = "corr_coeff")
  tc_p <- melt(tc_p, varnames = c("trait_1", "trait_2"), value.name = "p_value")

  # Merge into one data frame
  tc <- merge(tc_c, tc_p, by = c("trait_1", "trait_2"))
  
  ## Add species:
  tc$species <- species
  
  ## Return:
  return(tc)

} ## End of function.

## Create list to store outputs:
cls <- list() 

## Run function for all species:
for (p in unique(dp$species)) {  cls[[p]] <- run_co(p) }

## Collapse list into dataframe:
cor_df <- do.call(rbind, cls)

## Keep only the lower triangle (including diagonal):
cor_df <- cor_df %>% filter(as.numeric(trait_1) <= as.numeric(trait_2))

## Define factors for neat plots:
cor_df$trait_1 <- factor(x = cor_df$trait_1, levels = traits_selected, labels = trait_labels)
cor_df$trait_2 <- factor(x = cor_df$trait_2, levels = traits_selected, labels = trait_labels)

## Define species factor in desired facet order:
cor_df$species <- factor(
  x = cor_df$species,
  levels = c("muralis", "empty", "filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"),
  labels = c("P. muralis", "", "P. filfolensis", "P. pityusensis", "P. vaucheri", "P. liolepis", "P. cretensis", "P. gaigeae"))

## Plot heat maps of trait correlation coefficients:
f2b <- ggplot(cor_df, aes(x = trait_1, y = trait_2)) +
  
  ## Correlation coefficients:
  geom_tile(aes(fill = corr_coeff), colour = "black", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdBu", direction = 1, type = "div", na.value = "gray80", limits = c(-1, 1)) +
  
  ## Add asterisks for significance:
  geom_text(data = subset(cor_df, p_value <= 0.05), aes(label = "*"), color = "black", size = 4, nudge_y = -0.15) +

  ## Plot appearance:
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev) +
  labs(fill = "Correlation coefficients") +
  facet_wrap(facets = "species", ncol = 2, drop = FALSE) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    strip.text = element_text(face = "italic", size = 12, colour = "black"),
    strip.background = element_rect(fill = "gray90", colour = "transparent"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.775, 0.915),
    panel.spacing = unit(0.75, "lines"))
    
## Check:
f2b

## 5. Combine plots to compose Fig. 2 ----

## Combine:
f2 <- (f2a | f2b) +

  ## Relative widths of subplots:
  plot_layout(widths = c(1, 1)) +
  
  ## Subplot tags:
  plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag.position = c(0, 0.985),
      plot.tag = element_text(size = 16, face = "bold"))

## Check: 
f2

## Save:
ggsave(plot = f2, filename = paste0(path, "figures/Fig2_forest.jpg"), width = 32, height = 29, units = "cm")
ggsave(plot = f2, filename = paste0(path, "figures/Fig2_forest.pdf"), width = 32, height = 29, units = "cm", device = cairo_pdf)

## 6. Run trajectory analyses ----

## Status:
print("Now estimating trajectories!")

## Remove missing data that might have lingered:
nm <- na.omit(dp)

## After removing individuals with NA for any trait, how many individuals from each species in dataset?
table(nm$species)

## ID variables:
ID_var <- c("unique_identifier", "species", "species_abbpop", "group")

## Remove ID variables:
mx <- nm[setdiff(names(nm), ID_var)]

## Keep trait subset:
traits_selected
mx <- mx[traits_selected]

## Traits as numeric:
mx <- sapply(mx, as.numeric)

## To dataframe:
mx <- as.data.frame(mx)

## Row names:
row.names(mx) <- nm$unique_identifier

## Save:
#write.csv(x = mx, row.names = TRUE, file = "data/trajectory_analysis_trait_values.csv")

## Scale matrix:
mx <- scale(mx)

## Setting up data for geomorph:
rrppd <- rrpp.data.frame(Y = as.matrix(mx),
                         species = as.factor(nm$species),
                         group = as.factor(nm$group))

## Fit RRPP model:
rrppa <- lm.rrpp(data = rrppd, iter = 1000, f1 = Y ~ group + species + group * species)

## Model summary:
summary(rrppa, angle.type = "deg")

## ANOVA:
anova.lm.rrpp(rrppa)
  
## Trajectory analysis:
ta <- trajectory.analysis(fit = rrppa, pca = TRUE, groups = nm$species, traj.pts = nm$group)
  
## Summary of TA:
summary(ta)

## Check variable loadings:
ta$pca$rotation

## Percentages of PCs:
summary(ta$pca)

## Correlations (angles) between trajectories:
tc <- summary(ta, attribute = "TC", angle.type = "deg", show.trajectories = TRUE)

## Magnitude difference (absolute difference between path distances):
md <- summary(ta, attribute = "MD", show.trajectories = TRUE)
  
## Combine trajectory parameters:  
md_df <- md$summary.table
tc_df <- tc$summary.table
md_df$comparison <- row.names(md_df)
tc_df$comparison <- row.names(tc_df)
ta_df <- merge(md_df, tc_df, by = "comparison")

## Edit columns:
names(ta_df) <- c("comparison", 
                  "distance", "distance_UCL95", "distance_Z", "distance_p", "r",
                  "angle", "angle_UCL95", "angle_Z", "angle_p")
ta_df <- ta_df %>% separate(comparison, into = c("species_1", "species_2"), sep = ":")

## Or, use raw magnitude differences (preserving directions, i.e., signs):
## Extract path distances for each species:
rd <- ta$PD$obs
rd <- outer(rd, rd, FUN = "-") ## Subtract magnitudes.
rd <- melt(rd, varnames = c("species_1", "species_2"), value.name = "rd") ## Long format.
rd <- rd[rd$species_1 != rd$species_2, ] ## Remove diagonal.
rd <- rd[as.integer(rd$species_1) < as.integer(rd$species_2), ] ## Upper triangle.

## Merge:
ta_df <- merge(ta_df, rd, by = c("species_1", "species_2"))
View(ta_df)

## 7. Plot trajectories ----

## Get data to plot:
tp <- plot(ta)

## Subset and organize to plot:
## PCA values:
pd <- data.frame(
  unique_identifier = nm$unique_identifier,
  species = rrppd$species,
  group = rrppd$group,
  PC1 = tp$pc.points[, 1],
  PC2 = tp$pc.points[, 2])

## List to store trajectory points and lines:
pt_p <- list()
pt_l <- list()

## Loop over species:
for (s in names(tp$trajectories)) {
  pt_p[[s]] <- data.frame(species = s,
                          group = c("Brown", "Green"),
                          PC1 = tp$trajectories[[s]][, 1],
                          PC2 = tp$trajectories[[s]][, 2])
  pt_l[[s]] <- data.frame(species = s,
                          x1 = tp$trajectories[[s]][1, ][1], 
                          y1 = tp$trajectories[[s]][1, ][2], 
                          x2 = tp$trajectories[[s]][2, ][1], 
                          y2 = tp$trajectories[[s]][2, ][2]) 
} ## End of loop (s).
  
## List to dataframe:
pt_p <- do.call(rbind, pt_p)
pt_l <- do.call(rbind, pt_l)

## Define factors for ordered plots:
pd$species <- factor(x = pd$species, levels = rev(taxa), labels = rev(taxon_labels))
pt_p$species <- factor(x = pt_p$species, levels = rev(taxa), labels = rev(taxon_labels))
pt_l$species <- factor(x = pt_l$species, levels = rev(taxa), labels = rev(taxon_labels))

## Plot PCAs:
f3a <- ggplot() +
    
  ## Points for all individuals:    
  geom_point(data = pd, alpha = 0.9, color = "gray30",
             size = 2, shape = 21, show.legend = FALSE,
             aes(y = PC2, x = PC1, fill = group)) +
  
  ## Points for trajectory centroids:
  geom_point(data = pt_p, alpha = 1, color = "black",
             size = 4, shape = 23, show.legend = FALSE,
             aes(y = PC2, x = PC1, fill = group)) +
  
  ## Lines for trajectories:
  geom_segment(data = pt_l, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black") +
  
  ## Facets:
  facet_wrap(~species, ncol = 3) +
  
  ## Individual labels:
  #geom_text(aes(label = unique_identifier), size = 2) +
  
  ## Appearance:
  scale_fill_manual(values = c("#9B5300", "#31B565")) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, colour = "black", face = "italic"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = 12))

## Check:
f3a

## Save:
ggsave(plot = f3a, width = 25, height = 25, units = "cm", filename = "plots/Fig3A.jpg")
ggsave(plot = f3a, width = 25, height = 25, units = "cm", filename = "plots/Fig3A.pdf", device = cairo_pdf)

## 8. Plot trajectory parameters ----

## Melt:
mta <- melt(ta_df, id.vars = c("species_1", "species_2"))

## Flip to get the upper triangle:
mta_flip <- mta %>%
  rename(sp1 = species_1, sp2 = species_2) %>%
  mutate(species_1 = sp2, species_2 = sp1) %>%
  select(species_1, species_2, variable, value)

## Combine both to get symmetric matrix:
mta_sym <- bind_rows(mta, mta_flip)

## Add diagonal (self comparisons):
all_species <- sort(unique(c(mta$species_1, mta$species_2)))
all_vars <- unique(mta$variable)
mta_diag <- expand.grid(species_1 = all_species,
                        species_2 = all_species,
                        variable = all_vars)

## Diagonals to NA or 0:
NA_val <- 0
#NA_val <- NA
mta_diag <- mta_diag %>% filter(species_1 == species_2) %>% mutate(value = NA_val)

## For self-comparisons, p should be 1 (i.e., no difference):
mta_diag$value[mta_diag$variable == "angle_p"] <- 1
mta_diag$value[mta_diag$variable == "distance_p"] <- 1

## Combine everything and remove duplicates:
mta_full <- bind_rows(mta_sym, mta_diag) %>%
  distinct(species_1, species_2, variable, .keep_all = TRUE)

## Apply factor ordering for plotting:
mta_full$species_1 <- factor(x = mta_full$species_1, levels = rev(taxa), labels = rev(taxon_labels))
mta_full$species_2 <- factor(x = mta_full$species_2, levels = rev(taxa), labels = rev(taxon_labels))

## Back to wide format for export:
mw <- pivot_wider(data = mta_full,
                  names_from = variable,
                  values_from = value,
                  id_cols = c(species_1, species_2))
write.csv(x = mw, file = "data/trajectory_analysis_angles_distances.csv", row.names = FALSE)

## For plotting, keep only the lower triangle (including diagonal):
mta_full_lower <- mta_full %>% filter(as.numeric(species_1) <= as.numeric(species_2))

## Check p-values:
#p_df <- mta_full_lower[mta_full_lower$variable == paste0("angle", "_p"), ]
#p_df <- mta_full_lower[mta_full_lower$variable == paste0("distance", "_p"), ]

## Function to plot pairwise trajectory attributes:
plot.ta <- function(var, var_max, legend_label, file_name) {
  
  ## Testing:
  #file_name <- "Fig3B"
  #file_name <- "Fig3_raw_distances"
  #var <- "rd"
  #var <- "distance"
  #var_max <- 2.5
  #legend_label <- "Raw distances"
  
  ## Data for plot:
  plot_df <- mta_full_lower[mta_full_lower$variable == var, ]
  
  ## P-value data:
  p_df <- mta_full_lower[mta_full_lower$variable == paste0(var, "_p"), ]
  
  ## Plot:
  pta <- ggplot(plot_df, aes(x = species_1, y = species_2)) +
  
  ## Plot tiles:
  geom_tile(aes(fill = value), colour = "grey70", linewidth = 0.2) +
  scale_fill_paletteer_c(limits = c(0, var_max), "grDevices::Blue-Yellow", direction = -1) +
  #scale_fill_distiller(palette = "RdBu", direction = 1, type = "div", na.value = "gray80", limits = c(-2.2, 2.2)) +                    
     
  ## Add asterisks for significance:
  geom_text(data = subset(p_df, value <= 0.05), 
            aes(label = "*"), color = "black", size = 4, nudge_y = -0.075) +
      
  ## Plot appearance:
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev) +
  labs(fill = legend_label) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black", face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, colour = "black", face = "italic"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.90, 0.75))
  
  ## Check:
  pta
  
  ## Save:
  #ggsave(plot = pta, width = 14, height = 12, units = "cm", filename = paste0("plots/", file_name, ".jpg"))
  
  return(pta)
  
} ## End of function.

## Plot:
f3b <- plot.ta(file_name = "Fig3B", var = "angle", var_max = 95, legend_label = "Angles")
f3c <- plot.ta(file_name = "Fig3C", var = "distance", var_max = 2.5, legend_label = "\u0394 Paths")

## 9. Combine plots to compose Fig. 3 ----

## Combine:
f3 <- f3a + ( f3b / f3c ) +
  
  ## Relative widths of subplots:
  plot_layout(widths = c(3, 1)) +
  
  ## Subplot tags:
  plot_annotation(tag_levels = 'A') &
                  theme(plot.tag.position = c(0, 1),
                        plot.tag = element_text(size = 16, face = "bold"),
                        plot.margin = margin(t = 10, r = 10, b = 5, l = 10))

## Check:
f3

## Save:
#ggsave(plot = f3, filename = paste0(path, "figures/Fig3.jpg"), width = 34, height = 20, units = "cm")
#ggsave(plot = f3, filename = paste0(path, "figures/Fig3.pdf"), width = 34, height = 20, units = "cm", device = cairo_pdf)

## End of script.
