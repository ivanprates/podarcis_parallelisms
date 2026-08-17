## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To estimate candidate gene composition dissimilarity across Podarcis species.
## To test whether dissimilarity predicts phenotypic trajectory attributes (angles, path distances).
## To compose manuscript figures.

## 1. Setting up ----

## Packages:
library(ape)
library(paletteer)
library(patchwork)
library(reshape2)
library(tidyverse)
library(vegan)

## Path to working directory: Ivan's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## Create folders to store outputs:
dir.create("gen_vs_morph_analyses/data/")
dir.create("gen_vs_morph_analyses/plots/")

## 2. Compile outlier gene lists ----

## Cutoff:
cutoff <- 0.9925

## Format:
cutoff <- formatC(cutoff, format = "f", digits = 4)

## Read shared outlier data:
genes <- read.csv(file = paste0("outliers/data_outliers/outlier_sharing_", cutoff, ".csv"), header = TRUE)

## Subset:
taxa <- sort(c("filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"))
g_pr <- genes[taxa]

## Edit labels:
names(g_pr) <- paste("P.", names(g_pr))
taxa <- names(g_pr)

## How many outlier genes per species?
apply(X = g_pr, MARGIN = 2, FUN = sum)

## 3. Compute candidate gene dissimilarity between species ----

## Matrix of the number of shared outlier genes among species pairs:
mat_shared <- t(g_pr) %*% as.matrix(g_pr)

## Estimate compositional dissimilarity:
## Sorensen, Jaccard:
## 1: No shared elements (completely dissimilar).
## 0: Identical sets (complete similarity).

## Sorensen:
mat_Sor <- vegdist(t(g_pr), method = "bray", binary = TRUE)
mat_Sor <- as.matrix(mat_Sor)

## Jaccard:
mat_Jac <- vegdist(t(g_pr), method = "jaccard")
mat_Jac <- as.matrix(mat_Jac)

## Dealing with diagonals (self-comparisons):
## Replace diagonals for NA if we don't want to including self-comparisons in the analyses:
diag(mat_shared) <- NA
diag(mat_Sor) <- NA
diag(mat_Jac) <- NA

## Long format for each index:
shared <- as.data.frame(mat_shared)
shared$species_1 <- rownames(shared)
shared <- pivot_longer(data = shared, cols = -species_1, names_to = "species_2", values_to = "index_value")
shared$index <- "N_shared_genes"

Sor <- as.data.frame(mat_Sor)
Sor$species_1 <- rownames(Sor)
Sor <- pivot_longer(data = Sor, cols = -species_1, names_to = "species_2", values_to = "index_value")
Sor$index <- "Sorensen"

Jac <- as.data.frame(mat_Jac)
Jac$species_1 <- rownames(Jac)
Jac <- pivot_longer(data = Jac, cols = -species_1, names_to = "species_2", values_to = "index_value")
Jac$index <- "Jaccard"

## Combine indexes:
shared <- rbind(shared, Sor, Jac)

## Save:
write.csv(shared, file = "gen_vs_morph_analyses/data/outlier_similarity_indexes.csv", row.names = FALSE)

## Plot heat maps of pairwise dissimilarity:
plot.hm <- function(var) {

  ## Testing:
  #var <- "Jaccard"
  
  ## Subset to focal index:
  shared_p <- shared
  shared_p <- shared_p[shared_p$index == var, ]
  
  ## Round digits:
  shared_p$value <- round(shared_p$index_value, 3)
  
  ## Define factors for ordering in plots:
  shared_p$species_1 <- factor(x = shared_p$species_1, levels = taxa)
  shared_p$species_2 <- factor(x = shared_p$species_2, levels = taxa)
  
  ## Keep only the lower triangle, including diagonal:
  #shared_p <- shared_p[as.numeric(shared_p$species_1) <= as.numeric(shared_p$species_2), ]

  ## Keep only the lower triangle, excluding diagonal:
  shared_p <- shared_p[as.numeric(shared_p$species_1) < as.numeric(shared_p$species_2), ]

  ## Plot:
  hm <- ggplot(shared_p, aes(x = species_1, y = species_2)) +
  
    ## Heatmap layer:
    geom_tile(aes(fill = index_value), colour = "grey70", linewidth = 0.2) +
    
    ## Axes:
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), limits = rev) +
    
    ## Colors:
    scale_fill_paletteer_c("grDevices::Blue-Yellow", direction = -1) +
    
    ## Labels:
    labs(fill = var) +
    
    ## Add values in each cell:
    #geom_text(aes(label = value), color = "black", size = 5) + 
    
    ## Appearance:
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(family = "Arial", color = "black", size = 12, face = "italic", angle = 45, hjust = 1),
      axis.text.y = element_text(family = "Arial", color = "black", size = 12, face = "italic"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "inside",
      legend.position.inside = c(0.90, 0.75))

  ## Save:
  ggsave(plot = hm, width = 18, height = 15, units = "cm",
         filename = paste0("gen_vs_morph_analyses/plots/pairwise_", var,".jpg"))
    
  ## Return:
  return(hm)

} ## Close function.

## Run to plot:
plot.hm("N_shared_genes")
plot.hm("Sorensen")
a <- plot.hm("Jaccard") ; a

## 4. Prepare trajectory data ----

## Read trajectory analysis results: Angles and path distances between species:
t <- read.csv(file = "morphological_analyses/data/trajectory_analysis_angles_distances.csv", header = TRUE)

## Subset columns:
t <- t[c("species_1", "species_2", "distance", "angle")]

## Column names:
names(t) <- c("species_1", "species_2", "Trajectory length difference", "Trajectory angle")

## Exclude comparisons to self (diagonals):
t <- t[t$species_1 != t$species_2, ]

## Prepare symmetric matrices:
## For the genetic results, these won't include muralis:
mat_angle <- acast(t, species_1 ~ species_2, value.var = "Trajectory angle")
mat_distance <- acast(t, species_1 ~ species_2, value.var = "Trajectory length difference")

## Create versions without muralis (since no genetic data from P. muralis was used):
mat_angle_nm <- mat_angle[!rownames(mat_angle) %in% "P. muralis", !colnames(mat_angle) %in% "P. muralis"]
mat_distance_nm <- mat_distance[!rownames(mat_distance) %in% "P. muralis", !colnames(mat_distance) %in% "P. muralis"]

## 5. Estimate phylogenetic distances ----

## Read tree:
tree <- read.tree("gen_vs_morph_analyses/tree/Yang_et_al_2021_wall_lizards_time-tree.nex")

## Ladderize:
tree <- ladderize(tree)

## Extract focal taxa:
ptaxa <- c("cretensis", "filfolensis", "gaigeae", "liolepis", "muralis_IT", "pityusensis", "vaucheri_SS")
tree <- keep.tip(tree, ptaxa)

## Rename tips to match the other matrices:
tree$tip.label <- gsub(x = tree$tip.label, pattern = "_.+", replacement = "")
tree$tip.label <- paste("P.", tree$tip.label)

## Check:
plot(tree)

## Cophenetic (patristic) distances:
phylo_dist <- cophenetic.phylo(tree)

## Reoder matrix to match the other matrices:
ptaxa <- sort(unique(t$species_1))
phylo_dist <- phylo_dist[ptaxa, ptaxa]

## Multiply to get phylo distances in millions of years:
phylo_dist <- phylo_dist * 100

## Save:
write.csv(phylo_dist, file = "gen_vs_morph_analyses/data/phylogenetic_distances_Yang_et_al.csv", row.names = FALSE)

## Create version without muralis (since no genetic data from muralis was used):
phylo_dist_nm <- phylo_dist[!rownames(phylo_dist) %in% "P. muralis", !colnames(phylo_dist) %in% "P. muralis"]

## Convert to long format:
p <- as.data.frame(as.table(phylo_dist), stringsAsFactors = FALSE)
names(p) <- c("species_1", "species_2", "Phylogenetic distance")

## Exclude comparisons to self (diagonals):
p <- p[p$species_1 != p$species_2, ]

## Merge phylo dist info to trajectory dist info:
tp <- merge(t, p, by = c("species_1", "species_2"))

## 6. Mantel tests ----

## Check species order in matrices:
rownames(mat_angle)
colnames(mat_angle)
rownames(mat_distance)
colnames(mat_distance)
rownames(phylo_dist)
colnames(phylo_dist)

rownames(mat_angle_nm)
colnames(mat_angle_nm)
rownames(mat_distance_nm)
colnames(mat_distance_nm)
rownames(phylo_dist_nm)
colnames(phylo_dist_nm)

## Run Mantel: Trajectories vs. candidate gene dissimilarity:
mantel(xdis = mat_Jac, ydis = mat_angle_nm, method = "pearson", permutations = 10000)
mantel(xdis = mat_Jac, ydis = mat_distance_nm, method = "pearson", permutations = 10000)

## Run Mantel: Phylo distances vs. trajectories and candidate gene dissimilarity:
mantel(xdis = mat_Jac, ydis = phylo_dist_nm, method = "pearson", permutations = 10000)
mantel(xdis = mat_angle, ydis = phylo_dist, method = "pearson", permutations = 10000)
mantel(xdis = mat_distance, ydis = phylo_dist, method = "pearson", permutations = 10000)

## 7. Bivariate plots ----

## Dissimilarity metrics to wide format:
s <- shared %>% pivot_wider(
  names_from = index, 
  values_from = index_value)

## Combine all metrics:
m <- merge(s, tp, by = c("species_1", "species_2"))

## Function:
plot.bi <- function(xvar, yvar) {

  ## Testing:
  #xvar <- "Jaccard"
  #yvar <- "Trajectory angle"
  #yvar <- "Jaccard"
  #xvar <- "Phylogenetic distance"
  
  ## Plot:
  bp <- ggplot(data = m, aes(y = .data[[yvar]], x = .data[[xvar]])) +
    
    ## Line:
    geom_smooth(method = "lm", se = FALSE, colour = "gray50") +
    
    ## Points:
    geom_point(size = 1.5, show.legend = FALSE) +
    
    ## Appearance:
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(family = "Arial", color = "black", size = 14),
      axis.text.x = element_text(family = "Arial", color = "black", size = 12),
      axis.text.y = element_text(family = "Arial", color = "black", size = 12),
      plot.margin = margin(l = 20))
  
  ## Return:
  return(bp)

} ## End of function.

## Select two focal variables and plot biplots:
b <- plot.bi(xvar = "Jaccard", yvar = "Trajectory angle") ; b
c <- plot.bi(xvar = "Jaccard", yvar = "Trajectory length difference") ; c
d <- plot.bi(xvar = "Phylogenetic distance", yvar = "Jaccard") ; d
e <- plot.bi(xvar = "Phylogenetic distance", yvar = "Trajectory angle") ; e
f <- plot.bi(xvar = "Phylogenetic distance", yvar = "Trajectory length difference") ; f

## 8. Compose figure (Fig. 5) ----

## Combine:
abc <- wrap_plots(a, b, c, ncol = 3, widths = c(1, 1, 1))
def <- wrap_plots(d, e, f, ncol = 3, widths = c(1, 1, 1))

## Combine all:
f5 <- wrap_plots(abc, def, nrow = 2) +  
  
  ## Tags:
  plot_annotation(tag_levels = 'A') &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 16, face = "bold"))


## Check:
f5

## Save:
ggsave(plot = f5, filename = paste0(path, "figures/FigS3_2026-08.jpg"), width = 36, height = 22, units = "cm")
ggsave(plot = f5, filename = paste0(path, "figures/FigS3_2026-08.pdf"), width = 36, height = 22, units = "cm", device = cairo_pdf)

## End of script.
