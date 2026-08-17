## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, March 2026.
## The goals of this R script are:
## To compose Fig. 1.

## Packages:
library(ape)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggtree)
library(patchwork)

## Path to working directory: IP's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## 1. Plot map ----

## Sampling sites:
s <- read.csv(file = "sample_information/sampling_sites.csv", header = TRUE)

## World:
world <- ne_countries(scale = "medium", returnclass = "sf")

## Extent of the Mediterranean:
med_bbox <- st_bbox(c(
  xmin = -10,
  xmax = 28,
  ymin = 26,
  ymax = 44), 
  crs = st_crs(world))

## Crop:
med <- st_crop(world, med_bbox)

## Plot map:
mp <- ggplot() +
  geom_rect(aes(xmin = -10, xmax = 28, ymin = 26, ymax = 44), fill = "azure") +
  geom_sf(data = med, fill = "grey85", color = "gray40", linewidth = 0.45) +
  coord_sf(crs = "EPSG:4326", xlim = c(-10, 28), ylim = c(30, 44), expand = FALSE) +
  geom_point(data = s, aes(x = lon, y = lat, fill = GB), color = "black", size = 5, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = c("#9B5300", "#31B565")) +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 1.5, colour = "gray30"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

mp## 2. Plot phylogeny ----

## Read tree:
tree <- read.tree("gen_vs_morph_analyses/tree/Yang_et_al_2021_wall_lizards_time-tree.nex")

## Ladderize:
tree <- ladderize(tree)

## Taxa:
sort(tree$tip.label)

## Remove redundant taxa:
rtaxa <- c("hispanicus_Gal", "muralis_IT", "muralis_SA", "muralis_SP", "OU2_AB", "peloponnesiaca_W", "siculus_S", "tiliguerta_S", "vaucheri_SS")
tree <- drop.tip(tree, rtaxa)

## Rename tips:
tree$tip.label <- gsub(x = tree$tip.label, pattern = "_.+", replacement = "")
tree$tip.label <- paste("P.", tree$tip.label)

## Taxa to bold:
btaxa <- c("P. pityusensis", "P. gaigeae", "P. cretensis", "P. filfolensis", "P. muralis", "P. vaucheri", "P. liolepis")

## Plot tree:
tp <- ggtree(tree, ladderize = FALSE) +
  geom_tiplab(family = "Arial", size = 6,
              aes(fontface = ifelse(label %in% btaxa, "bold.italic", "italic")))

## Expand axis:
te <- tp + xlim(0, max(tp$data$x) + 0.2)
te

## 3. Combine plots ----

## Combine:
l <- mp / plot_spacer() + plot_layout(heights = c(5, 1))
f1 <- wrap_plots(l, te, ncol = 2, widths = c(5, 2))
f1

## Save:
ggsave(plot = f1, filename = paste0(path, "figures/Fig1.jpg"), width = 35, height = 20, units = "cm")
ggsave(plot = f1, filename = paste0(path, "figures/Fig1.pdf"), width = 35, height = 20, units = "cm", device = cairo_pdf)

## End of script.
