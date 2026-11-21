## R script: Comparative dorsal coloration analyses across 7 Podarcis species
## Expanding Javier Abalos' (JA) suggestion to the full empirical dataset
## Author: Ivan Prates et al.
## Date: November 2026

## 1. Setup and Libraries ----
library(ggplot2)
library(patchwork)

options(scipen = 999)

## Define paths:
base_dir <- "/home/ivan/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms"
work_dir <- file.path(base_dir, "drafts/3_version_incorporating_coauthor_feedback")
data_path <- file.path(base_dir, "morphological_analyses/data/Podarcis_greenbrowns_LizWiz_outputs_sample_info_2025-06-25.csv")
meta_path <- file.path(work_dir, "Table_S1_sample_information.csv")

setwd(work_dir)

## 2. Load and Filter Data ----
df_lizwiz <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1")
table_s1 <- read.csv(meta_path, header = TRUE, stringsAsFactors = FALSE)

## Filter to dorsal coloration slice (dcc_x1-col2 / dorsum_total_color) for the 210 sequenced study specimens:
dors <- df_lizwiz[df_lizwiz$slice_label == "dorsum_total_color" & df_lizwiz$unique_identifier %in% table_s1$unique_identifier, ]

## Species ordering matching manuscript:
species_levels <- c("muralis", "filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae")
species_labels <- c("P. muralis", "P. filfolensis", "P. pityusensis", "P. vaucheri", "P. liolepis", "P. cretensis", "P. gaigeae")

dors$species_factor <- factor(dors$species, levels = species_levels, labels = species_labels)
dors$group_factor <- factor(dors$group, levels = c("Brown", "Green"))

cat(sprintf("Loaded %d individuals across %d species.\n", nrow(dors), length(unique(dors$species))))

## 3. Compute Summary Statistics ----
calc_stats <- function(var_name) {
  res <- do.call(rbind, lapply(split(dors, list(dors$species_factor, dors$group_factor)), function(sub) {
    if (nrow(sub) == 0) return(NULL)
    val <- sub[[var_name]]
    data.frame(
      trait = var_name,
      species = sub$species_factor[1],
      group = sub$group_factor[1],
      n = length(val),
      mean = mean(val, na.rm = TRUE),
      sd = sd(val, na.rm = TRUE),
      median = median(val, na.rm = TRUE),
      min = min(val, na.rm = TRUE),
      max = max(val, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(res) <- NULL
  return(res)
}

stats_hue <- calc_stats("h.theta")
stats_sat <- calc_stats("r.achieved")
stats_lum <- calc_stats("lum_s")

summary_stats <- rbind(stats_hue, stats_sat, stats_lum)
write.csv(summary_stats, "dorsal_coloration_summary_stats_7_species.csv", row.names = FALSE)
cat("Summary statistics saved to dorsal_coloration_summary_stats_7_species.csv\n")

## 4. Pairwise Interspecific Hue Differences ----
means_br <- tapply(dors$h.theta[dors$group == "Brown"], dors$species[dors$group == "Brown"], mean)
means_gr <- tapply(dors$h.theta[dors$group == "Green"], dors$species[dors$group == "Green"], mean)

pw_list <- list()
idx <- 1
for (i in 1:(length(species_levels) - 1)) {
  for (j in (i + 1):length(species_levels)) {
    sp1 <- species_levels[i]
    sp2 <- species_levels[j]
    d_br <- abs(means_br[sp1] - means_br[sp2])
    d_gr <- abs(means_gr[sp1] - means_gr[sp2])
    
    pw_list[[idx]] <- data.frame(
      pair = paste(sp1, sp2, sep = " - "),
      sp1 = sp1,
      sp2 = sp2,
      diff_brown = d_br,
      diff_green = d_gr,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}
pw_df <- do.call(rbind, pw_list)
write.csv(pw_df, "dorsal_coloration_pairwise_contrasts_7_species.csv", row.names = FALSE)
cat("Pairwise contrasts saved to dorsal_coloration_pairwise_contrasts_7_species.csv\n")

## 5. Statistical Modeling ----
## Linear model testing interaction of species and phenotype on hue:
lm_hue <- lm(h.theta ~ species_factor * group_factor, data = dors)
aov_hue <- anova(lm_hue)
cat("\nANOVA for dorsal hue (h.theta):\n")
print(aov_hue)

## Variance test (Levene-type / F-test across all brown vs all green):
f_test_hue <- var.test(dors$h.theta[dors$group == "Green"], dors$h.theta[dors$group == "Brown"])
cat("\nF-test for equality of variance in hue (Green vs Brown):\n")
print(f_test_hue)

## 6. Plotting Fig. S2 ----
color_palette <- c("Brown" = "#9B5300", "Green" = "#31B565")

theme_ms <- theme_classic() +
  theme(
    text = element_text(family = "sans"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9, colour = "black", face = "italic"),
    axis.text.y = element_text(size = 9, colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )

## Panel A: Dorsal Hue (theta)
p_hue <- ggplot(dors, aes(x = species_factor, y = h.theta, fill = group_factor, color = group_factor)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, alpha = 0.25, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  labs(y = expression(bold("Dorsal hue (")*bold(theta)*bold(", radians)")), title = "Dorsal Hue") +
  theme_ms

## Panel B: Dorsal Saturation (r.achieved)
p_sat <- ggplot(dors, aes(x = species_factor, y = r.achieved, fill = group_factor, color = group_factor)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, alpha = 0.25, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  labs(y = expression(bold("Dorsal saturation (")*bolditalic(r)[achieved]*bold(")")), title = "Dorsal Saturation") +
  theme_ms

## Panel C: Dorsal Luminance (lum_s)
p_lum <- ggplot(dors, aes(x = species_factor, y = lum_s, fill = group_factor, color = group_factor)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, alpha = 0.25, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  labs(y = expression(bold("Dorsal luminance (")*bolditalic(L)[s]*bold(")")), title = "Dorsal Luminance") +
  theme_ms

## Panel D: Interspecific Pairwise Differences in Mean Hue
pw_long <- data.frame(
  phenotype = rep(c("Brown morphs", "Green morphs"), each = nrow(pw_df)),
  diff = c(pw_df$diff_brown, pw_df$diff_green)
)
pw_long$phenotype <- factor(pw_long$phenotype, levels = c("Brown morphs", "Green morphs"))

p_pw <- ggplot(pw_long, aes(x = phenotype, y = diff, fill = phenotype, color = phenotype)) +
  geom_boxplot(width = 0.45, alpha = 0.25, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.75) +
  scale_fill_manual(values = c("Brown morphs" = "#9B5300", "Green morphs" = "#31B565")) +
  scale_color_manual(values = c("Brown morphs" = "#9B5300", "Green morphs" = "#31B565")) +
  labs(y = expression(bold("Interspecific |")*bold(Delta)*bold(theta)*bold("| (radians)")), title = "Interspecific Hue Divergence") +
  theme_classic() +
  theme(
    text = element_text(family = "sans"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 9.5, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 9, colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )

## Combine with patchwork:
fig_s2 <- (p_hue | p_sat) / (p_lum | p_pw) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14, face = "bold"))

## Save Figure:
ggsave("Fig_S2_dorsal_coloration.pdf", plot = fig_s2, width = 24, height = 20, units = "cm")
ggsave("Fig_S2_dorsal_coloration.png", plot = fig_s2, width = 24, height = 20, units = "cm", dpi = 300)

cat("Figures saved to Fig_S2_dorsal_coloration.pdf and Fig_S2_dorsal_coloration.png\n")
