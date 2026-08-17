## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To employ simulations to whether the number of overlapping candidate genes 
## across Podarcis species is higher than expected by chance.

## Packages:
library(tidyverse)

## 1. Getting ready ----

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/")

## Create directories to save outputs:
dir.create("outlier_sharing_simulations/")
dir.create("outlier_sharing_simulations/plots")
dir.create("outlier_sharing_simulations/data")

## Define function: Extract gene IDs for a species:
get_species_genes <- function(df, sp_col) {
  unique(as.character(df$ID)[df[[sp_col]] == 1])
}

## Define function: simulated overlapping gene sets:
simulate_overlap <- function(n_outliers, all_genes) {
  sampled_genes <- lapply(n_outliers, function(n) sample(all_genes, n))
  sim_overlap <- Reduce(intersect, sampled_genes)
  length(sim_overlap)
}

## 3. Perform simulations under each cutoff ----

## Cutoff:
cutoff <- 0.9925

## Or, loop over cutoffs:
#for (cutoff in c(0.9900, 0.9925, 0.9950)) {

## Format:
cutoff <- formatC(cutoff, format = "f", digits = 4)

## Read data:
genes <- read.csv(file = paste0("outliers/data_outliers/outlier_sharing_", cutoff, ".csv"), header = TRUE)

## Remove column:
#genes <- genes[setdiff(names(genes), "Pmu_candidate")]

## Or rather, rename column:
names(genes)[names(genes) == "Pmu_candidate"] <- "muralis"

## Species considered:
species <- setdiff(names(genes), c("ID", "gene", "shared_by"))
species

## Genes in P. muralis:
table(genes$muralis)
genes[genes$muralis == 1, ]

## Read NCC genes from Feiner et al. 2024:
NCC <- read.csv(file = "baypass/annotations/FST_outliers_NCC_genes.csv")
NCC$gene <- toupper(NCC$gene)
sum(NCC$NCC_gene)

## NCC genes in outliers:
NCC$gene[NCC$gene %in% genes$gene]

## Limit gene set to all muralis genes:
genes <- genes[genes$muralis == 1, ]

## Check for species:
sum(genes$cretensis)
sum(genes$gaigeae)
sum(genes$filfolensis)
sum(genes$liolepis)
sum(genes$pityusensis)
sum(genes$vaucheri)

## Alternatively, limit gene set to NCC genes:
#genes <- genes[genes$gene %in% NCC$gene, ]

## All P. muralis genes:
Pmu_genes <- read.table(file = "baypass/annotations/GCF_004329235.1_PodMur_1.0_genomic_patched_candidate_region_genes.bed", header = FALSE)
names(Pmu_genes) <- c("CHROM", "start", "end", "gene", "ID")

## Keep only gene IDs:
Pmu_genes <- Pmu_genes$ID

## Dataframe to store results:
sim_df <- data.frame()

## Loop over k species:
for (k in 2:length(species)) {

  ## Testing:
  #k <- 2

  ## Status:
  #cat("Now simulating overlaps for sets of", k, "species!\n")
  
  ## Generate all combinations of k species:
  sp_comb <- combn(x = species, m = k, simplify = FALSE)
  
  ## Loop over species combinations:
  for (cb in sp_comb) {
  
    ## Testing:
    #cb <- sp_comb[[1]]
    
    ## Status:
    print(paste("Cutoff:", cutoff))
    print(cb)
    
    ## Number of outlier genes per species:
    n_outliers <- sapply(cb, function(sp) sum(genes[[sp]] == 1, na.rm = TRUE))
    
    ## Extract and combine genes across species in the combination:
    gene_ls <- lapply(cb, function(sp) get_species_genes(genes, sp))
    
    ## List empirical overlapping genes in the combination:
    emp_overlap <- Reduce(intersect, gene_ls)
    
    ## Number of empirical overlapping genes:  
    n_emp <- length(emp_overlap)
    
    ## If the number of empirical overlapping genes is not zero:
    if (n_emp >= 1) {
      
    ## Simulations:
    set.seed(123)
    n_sim <- 10000
    sim_overlaps <- replicate(n = n_sim, expr = simulate_overlap(n_outliers = n_outliers, all_genes = Pmu_genes))
    
    ## One-sided p-value (excess overlaps):
    p_value <- mean(x = sim_overlaps >= n_emp)
    
    ## Save results:
    results <- data.frame(
      species_combination = paste(cb, collapse = "-"),
      No_species_in_set = k,
      No_genes_shared_empirical = n_emp,
      No_genes_shared_null_median = median(sim_overlaps),
      p_value = p_value)
    sim_df <- rbind(sim_df, results)
    
    ## Plot null distributions:
    png(filename = paste0("outlier_sharing_simulations/plots/sims_", cutoff, "_", cutoff, "_", paste(cb, collapse = "-"), ".jpg"),
        width = 20, height = 15, units = "cm", res = 150)
    hist(sim_overlaps,
         main = paste0("Null overlap distribution among ", k, " species:"),
         xlab = "Number of overlapping genes",
         col = "lightblue", border = "black",
         xlim = range(c(sim_overlaps, n_emp)))
    mtext(paste(cb, collapse = ", "))
    abline(v = n_emp, col = "red", lwd = 2)
    legend(x = "right", bty = "n", legend = paste("Obs. overlaps:", n_emp, 
                                                  "\nMedian of null:", median(sim_overlaps), 
                                                  "\np-value:", signif(p_value, 3)))
    dev.off()
    
    } ## Close if statement (n_emp).
  
  } ## Close loop (cb).
  
} ## Close loop (k).

## Write results to file:
#write.csv(sim_df, paste0("outlier_sharing_simulations/data/sims_", cutoff, ".csv"), row.names = FALSE)

## Limit to comparisons with muralis:
mur_df <- sim_df[sim_df$No_species_in_set == 2, ]
mur_df <- separate(data = mur_df, col = species_combination, into = c("sp1", "sp2"))
mur_df <- mur_df[mur_df$sp1 == "muralis" | mur_df$sp2 == "muralis", ]

## Check:
View(mur_df)

#} ## Close loop (cutoff).

## End of script.
