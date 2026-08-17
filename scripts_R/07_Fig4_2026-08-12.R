## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To compose manuscript figure 4.

## 1. Getting ready ----

## Packages:
library(ComplexUpset)
library(cowplot)
library(ggplot2)
library(grid)
library(jpeg)
library(patchwork)
library(tidyverse)

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Path to working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## 2. Data ----

## Taxa:
taxa <- sort(c("filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"))

## Cutoff:
cutoff <- 0.9925

## Or, loop over combinations:
#for (cutoff in c(0.9900, 0.9925, 0.9950)) {

## Format:
cutoff <- formatC(cutoff, format = "f", digits = 4)

## Data:
g_sh <- read.csv(file = paste0("outliers/data_outliers/outlier_sharing_", cutoff, ".csv"), header = TRUE)

## Exclude columns:
g_pr <- g_sh[taxa]

## Rename columns:
names(g_pr) <- paste("P.", names(g_pr), sep = " ")

## 3. Outlier genes shared by N species ----

## Count:
gene_counts <- rowSums(g_pr)

## Tabulate exact counts:
gene_counts <- as.data.frame(table(n_species = gene_counts))

## Barplot:
f <- ggplot(gene_counts, aes(x = n_species, y = Freq)) +
  geom_col(fill = "gray80", color = "black") +
  geom_text(aes(label = Freq), vjust = -0.5, family = "Arial", size = 4, colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "No. species sharing outlier", y = "No. outliers") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.line = element_line(),
    axis.ticks = element_line(),
    panel.grid = element_blank())

## Check:
f

## 4. No. outlier genes per species ----
species_counts <- data.frame(species = colnames(g_pr), n_outlier_genes = colSums(g_pr))

## Plot:
g <- ggplot(species_counts, aes(x = reorder(species, n_outlier_genes), y = n_outlier_genes)) +
  geom_col(fill = "gray30") +
  coord_flip() +
  #scale_y_reverse() +
  labs(x = NULL, y = "No. outlier genes per species") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(face = "italic", size = 14, colour = "black", hjust = 0.5),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black"),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.x = element_line())

## Check:
g

## 5. Upset plots of outlier gene sharing ----

## Upset plot:
de <- upset(
  data = g_pr,
  intersect = colnames(g_pr),
  mode = "inclusive_intersection",
  #stripes = "white",
  name = "Species set",
  #width_ratio = 0.15,
  height_ratio = 0.8,
  
  ## Editing set size histogram:
  set_sizes = FALSE, 
  #set_sizes = (upset_set_size() +
  #             labs(y = "No. outlier genes per species") +
  #             theme(
  #               text = element_text(family = "Arial"),
  #               axis.line.x = element_line(),
  #               axis.ticks.x = element_line(),
  #               panel.grid.major.x = element_blank(),
  #               panel.grid.minor.x = element_blank())),
  
  ## Editing intersection histogram:
  base_annotations = list(
        "No. genes in set" = intersection_size(counts = TRUE) +
        theme(
          text = element_text(family = "Arial"),
          #axis.text.y = element_text(size = 10, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()))) +
  
  ## Editing sets plot:
  theme(
    axis.title.x = element_text(size = 14, family = "Arial", colour = "black"),
    #axis.text.y = element_text(size = 14, face = "italic", family = "Arial", colour = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()) 

## Check:
de

## 6. Row indicating significance:

## Upset data:
ud <- upset_data(g_pr, intersect = colnames(g_pr), mode = "inclusive_intersection")

## Intersection names in plotting order:
intersections <- ud$sorted$intersections

## Significance data:
sd <- read.csv(file = paste0("outlier_sharing_simulations/data/sims_", cutoff, ".csv"), header = TRUE)

## Simplify:
sd <- sd[c("species_combination", "p_value")]

## Add rows for single species:
ss <- data.frame(species_combination = taxa, p_value = 1)
sd <- rbind(sd, ss)

## Split combinations into species lists:
sd$species_list <- strsplit(sd$species_combination, "-")

## Keep only combinations composed entirely of taxa in the upset plot:
sd <- sd[sapply(sd$species_list, function(x) all(x %in% taxa)), ]

## Presence/absence matrix:
presence <- t(sapply(sd$species_list, function(x) taxa %in% x))
colnames(presence) <- taxa

## Build intersection IDs exactly as in ComplexUpset:
make_id <- function(x) { paste(paste0("P. ", taxa[x]), collapse = "-") }
sd$intersection_id <- apply(presence, 1, make_id)

## Intersections plotted:
sd <- sd[sd$intersection_id %in% intersections, ]

## Enforce plotting order:
sd$intersection_id <- factor(sd$intersection_id, levels = rev(intersections))

## Plot asterisks for significance:
p_sig <- ggplot(data = sd, 
          aes(x = intersection_id, y = 1, label = ifelse(p_value < 0.05, "*", ""))) +
  geom_text(size = 6, vjust = 0.5) +
  theme_void() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, 0, 0, 0),
        plot.background = element_blank())

## Check:
p_sig

## 4. Compose figure (Fig. 4) ----

## Add tags to each plot manually:
p_sig <- p_sig + labs(tag = "A")
f <- f + labs(tag = "B")
g <- g + labs(tag = "C")

## Compose:
fg <- wrap_plots(f, g, nrow = 2, heights = c(1.25, 1))
dex <- wrap_plots(p_sig, de, nrow = 2, heights = c(1, 10))
f4 <- wrap_plots(dex, fg, ncol = 2, widths = c(5, 1)) &
  theme(
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 16, face = "bold"))

## Save:
ggsave(plot = f4, width = 40, height = 13, units = "cm", dpi = 300, 
  filename = paste0(path, "figures/Fig4_c", cutoff, ".jpg"))
ggsave(plot = f4, width = 40, height = 13, units = "cm", device = cairo_pdf,
  filename = paste0(path, "figures/Fig4_c", cutoff, ".pdf"))

#} ## Close loop (cutoff).

## End of script.
