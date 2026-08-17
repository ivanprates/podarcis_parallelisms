## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To dissect a candidate region on chromosome 12.

## 1. Getting ready ----

## Packages:
library(patchwork)
library(tidyverse)

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Path to working directory: Ivan's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## Create dir:
dir.create("candidate_region/data")
dir.create("candidate_region/plots")

## Taxa:
taxa <- sort(c("Pfi", "Ppi", "Pva", "Pli", "Pcr", "Pga"))

## Candidate region:
region_start <- 21226312
region_end <- 21841019

## x limits for plots given region coordinates:
xlims <- c(region_start, region_end) / 1e6

## 2. Extract chrom 12 data ----

## Loop over taxa: 
# for (sp in taxa) {
#   
#   ## Testing:
#   #sp <- taxa[1]
#   
#   ## Status:
#   print(paste0("Now processing ", sp, "!"))
#   
#   ## Read XtX estimates alone for chrom 12:
#   merged <- read.table(file = paste0("baypass/MAF005/data_XtX/", sp, "_median_summary_pi_xtx.out"), header = TRUE)
#   
#   ## And then read SNP positional information:
#   SNP_ids <- read.table(file = paste0("baypass/MAF005/data_XtX/", sp, "_snp_ids.txt"), header = FALSE)
#   names(SNP_ids) <- c("CHROM", "SNP", "MRK")
#   
#   ## Merge both:
#   scan <- merge(SNP_ids, merged, by = "MRK")
#   
#   ## Clean initial objects:
#   rm(merged, SNP_ids)
#   
#   ## Focus on chrom 12:
#   scan <- scan[scan$CHROM == "NC_041323.1", ]
#   
#   ## Focus on presumed candidate region (614,707 bp long):
#   scan <- scan[scan$SNP >= region_start & scan$SNP <= region_end, ]
#   
#   ## Save it:
#   write.table(scan, file = paste0("candidate_region/data/", sp, "_median_summary_pi_xtx_SNP_ids_CHROM12_candidate_region.out"), row.names = FALSE)
# 
# } ## Close loop (sp).

## 3. Combine XtX data ----

## Window size for plots:
## This is for tractability, otherwise plotting might become impossible.
win_XtX <- 1 ## Not window, use all SNPs.

## List to store:
f_ls <- list()

## Loop over taxa:
for (sp in taxa) {
 
  ## Testing:
  #sp <- taxa[1]
  
  ## Read data:
  s_df <- read.table(paste0("candidate_region/data/", sp, "_median_summary_pi_xtx_SNP_ids_CHROM12_candidate_region.out"), header = TRUE)
  
  ## Define windows:
  windows <- s_df
  windows$window <- floor( ( windows$SNP - region_start ) / win_XtX ) + 1
  windows$win_start <- ( region_start + ( win_XtX * ( windows$window - 1 ) ) )
  windows$win_end <- windows$win_start + win_XtX - 1
  
  ## Estimate average XtX per window:
  windows <- windows %>%
    group_by(window, win_start, win_end) %>%
    summarise(XtX = mean(XtX, na.rm = TRUE), .groups = "drop")
  
  ## Mid position of windows:
  windows$win_mid <- (( windows$win_start + windows$win_end) / 2 )
  
  ## Divide position by 10⁶:
  windows$win_mid <- windows$win_mid/1000000

  ## Taxon:
  windows$taxon <- sp
  
  ## Store:
  f_ls[[sp]] <- windows

} ## Close loop (sp).

## Combine:
f_df <- do.call(rbind, f_ls)

## Ordering:
f_df$taxon <- factor(f_df$taxon, levels = taxa, 
                     labels = c("P. cretensis", "P. filfolensis", "P. gaigeae", 
                                "P. liolepis", "P. pityusensis", "P. vaucheri"))

## Now read cutoffs:
## List to store:
c_ls <- list()

## Loop over taxa:
for (sp in taxa) {
  #sp <- taxa[1]
  
  ## List to store:
  t_ls <- list()
  
  ## Loop over cutoffs:
  #for (cutoff in c(0.9900, 0.9925, 0.9950)) {
  
    ## Cutoff:
    cutoff <- 0.9925
    
    ## Format:
    cutoff <- formatC(cutoff, format = "f", digits = 4)
    
    ## Read, edits:
    c_df <- read.table(paste0("baypass/0326/data_XtX/", sp, "_mean_xtx_threshold_", cutoff, ".txt"), header = FALSE)
    names(c_df) <- "XtX_threshold"
    c_df$taxon <- sp
    c_df$Cutoff <- cutoff
    
    ## Store:
    t_ls[[cutoff]] <- c_df
  
  #} ## Close loop (c).

  ## Combine, store:
  c_ls[[sp]] <- do.call(rbind, t_ls)
  
} ## Close loop (sp).

## Combine, store:
c_df <- do.call(rbind, c_ls)

## labels:
c_df$taxon <- factor(c_df$taxon, levels = taxa, 
                     labels = c("P. cretensis", "P. filfolensis", "P. gaigeae", 
                                "P. liolepis", "P. pityusensis", "P. vaucheri"))

## 4. Plot XtX for candidate region ----

## Taxon sets to plot in pairs:
taxon_sets <- list(
  c("P. cretensis", "P. filfolensis"),
  c("P. gaigeae", "P. liolepis"),
  c("P. pityusensis", "P. vaucheri"))

## List to store plots:
taxon_XtX <- list()

## Loop over taxon sets:
for (i in seq_along(taxon_sets)) {
  
  ## Taxon set to plot:
  taxon_set <- taxon_sets[[i]]

  ## Plot:
  a <- ggplot(
    data = filter(f_df, taxon %in% taxon_set)) +

    ## Points:
    geom_point(
      mapping = aes(x = win_mid, y = XtX), 
      size = 0.5, color = "black") +

    ## Wrap:
    facet_wrap(~taxon, nrow = 1) +

    ## Cutoff lines:
    geom_hline(
      data = filter(c_df, taxon %in% taxon_set), 
      alpha = 0.5, show.legend = FALSE, colour = "blue",
      mapping = aes(yintercept = XtX_threshold)) +

    ## Labels:
    labs(
      x = "Position on chromosome 12 (x 10⁶)",
      y = "XtX statistic") +

    ## Appearance:
    scale_x_continuous(
      limits = xlims,
      expand = c(0, 0),
      breaks = c(21.25, 21.38, 21.52, 21.66, 21.8)) +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      text = element_text(family = "Arial", color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(colour = "black", linewidth = 0.5))

  ## Store in list:
  taxon_XtX[[i]] <- a

} ## End of loop (taxon_set).

## 5. Coverage for candidate region ----

## List to store:
d_ls <- list()

## Loop over taxa:
for (sp in taxa) {
  
  ## Testing:
  #sp <- taxa[1]
  
  ## Status:
  print(paste0("Now processing ", sp, "!"))
  
  ## Data:
  t_df <- read.table(paste0("candidate_region/depth/", sp, "_NC_041323_candidate_region.depth.txt"), header = FALSE)
  
  ## Names and color assignments per taxon:
  if (sp %in% c("Pfi", "Ppi", "Pva")) { 
    names(t_df) <- c("CHROM", "POS", paste0("G", 1:15), paste0("B", 1:15)) }
  if (sp %in% c("Pcr", "Pga", "Pli")) { 
    names(t_df) <- c("CHROM", "POS", paste0("B", 1:15), paste0("G", 1:15)) }
  
  ## Focus on presumed candidate region (614,707 bp long):
  t_df$POS <- as.numeric(t_df$POS)
  t_df <- t_df[t_df$POS >= region_start & t_df$POS <= region_end, ]

  ## Average depth across G samples:
  G_col <- grep(x = names(t_df), pattern = "G")
  t_df$avgDP_G <- rowMeans(t_df[, G_col], na.rm = TRUE)
  
  ## Average depth across B samples:
  B_col <- grep(x = names(t_df), pattern = "B")
  t_df$avgDP_B <- rowMeans(t_df[, B_col], na.rm = TRUE)
  
  # Cumulative depth across G samples:
  t_df$cumDP_G <- rowSums(t_df[, G_col], na.rm = TRUE)
  
  # Cumulative depth across B samples:
  t_df$cumDP_B <- rowSums(t_df[, B_col], na.rm = TRUE)
  
  ## Subset:
  t_df <- t_df[c("POS", "avgDP_G", "avgDP_B")]
  
  ## Add taxon, store:
  t_df$taxon <- sp
  d_ls[[sp]] <- t_df
  
} ## Close loop (sp).

## Combine, store:
d_df <- do.call(rbind, d_ls)

## Check:
head(d_df)

## Window size for plots:
## This is for tractability, otherwise plotting becomes impossible.
win_DP <- 100

## Define windows:
w_df <- d_df
w_df$window <- floor( ( w_df$POS - region_start ) / win_DP ) + 1
w_df$win_start <- ( region_start + ( win_DP * ( w_df$window - 1 ) ) )
w_df$win_end <- w_df$win_start + win_DP - 1

## Average depth per window:
w_df <- w_df %>%
  group_by(taxon, window, win_start, win_end) %>%
  reframe(
    mid = mean((win_start + win_end)/2),
    mean_avgDP_B = mean(avgDP_B, na.rm = TRUE),
    mean_avgDP_G = mean(avgDP_G, na.rm = TRUE))

## Divide position by 10⁶:
w_df$mid <- w_df$mid/1000000

## Ordering:
w_df$taxon <- factor(w_df$taxon, levels = taxa)

## Ordering:
w_df$taxon <- factor(w_df$taxon, levels = taxa, 
                     labels = c("P. cretensis", "P. filfolensis", "P. gaigeae", 
                                "P. liolepis", "P. pityusensis", "P. vaucheri"))

## Check:
head(w_df)

## List to store plots:
taxon_depth <- list()

## Loop over taxon sets:
for (i in seq_along(taxon_sets)) {
  
  ## Taxon set to plot:
  taxon_set <- taxon_sets[[i]]

  ## Plot:
  b <- ggplot(
    data = filter(w_df, taxon %in% taxon_set)) +
    
    ## Lines:
    geom_smooth(
      aes(x = mid, y = mean_avgDP_B), 
      method = "loess", span = 0.075, se = F, color = "#9B5300", linewidth = 1) +
    geom_smooth(
      aes(x = mid, y = mean_avgDP_G), 
      method = "loess", span = 0.075, se = F, color = "#31B565", linewidth = 1) +
    
    ## Facets:
    facet_wrap(~ taxon, nrow = 1) +
    
    ## Labels:
    labs(
      y = "Average read depth (x)") +
    
    ## Appearance:
    scale_x_continuous(limits = xlims, expand = c(0, 0)) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "italic", size = 14),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(family = "Arial", color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14),
      plot.title = element_text(color = "black", size = 16),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(colour = "black", linewidth = 0.5))

  ## Store in list:
  taxon_depth[[i]] <- b

} ## End of loop (taxon_set).

## 6. Bar with gene positions ----

## Gene positions:
genes <- read.csv("baypass/annotations/CandidateRegion_Genes_Location_Spa_exons.csv")
head(genes)

## Remove, following Feiner et al. 2023:
genes <- genes[genes$name != "ptchd3_X1", ]

## Capitalize:
genes$name <- str_to_title(genes$name)

## Edit:
genes$name <- gsub(x = genes$name, pattern = "_", replacement = "-")

## Combine exons into genes:
genes <- genes %>%
  group_by(name) %>%
  reframe(
    start = min(start),
    end = max(end),
    strand = first(strand))

## Adjust coordinates:
genes <- genes %>%
  mutate(
    xmin = start / 1e6,
    xmax = end / 1e6,
    strand_ymin = ifelse(strand == "+", 0, -0.1),
    strand_ymax = ifelse(strand == "+", 0.1, 0),
    label_y = ifelse(strand == "+", 0.2, -0.5))

## Plot gene bars:
g <- ggplot(genes) +

  ## Center line:
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray20") +

  ## Gene bars:
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = strand_ymin, ymax = strand_ymax),
    show.legend = FALSE, fill = "gray20", colour = "gray20") +

  ## Y-axis to fit gene labels:
  scale_y_continuous(limits = c(-0.7, 0.8), expand = c(0, 0)) +

  ## Gene labels:
  geom_text_repel(
    aes(x = (xmin + xmax) / 2, y = label_y, label = name),
    family = "Arial", fontface = "italic", color = "black", angle = 70, size = 4, 
    direction = "x", box.padding = 0.1, point.padding = 0, min.segment.length = Inf,
    hjust = 0) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = xlims, expand = c(0, 0)) +
  
  ## Looks:
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_blank(),
    plot.margin = margin(l = 3, r = 3))

## 7. Combine plots ----

## Loop over taxon sets:
for (i in 1:3) {
  
  ## Combine plots
  abg <- ( taxon_depth[[i]] / ( g + g ) / taxon_XtX[[i]] ) + plot_layout(heights = c(3, 3, 3))

  ## Save:
  ggsave(
    plot = abg,
    width = 30, height = 20, units = "cm", dpi = 300, 
    filename = paste0("figures/FigS2_", i, ".jpg"))
  
} ## Close loop (taxon_set)

## End of script.
