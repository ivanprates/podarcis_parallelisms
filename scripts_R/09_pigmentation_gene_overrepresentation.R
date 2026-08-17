## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To check whether pigmentation-related genes are overrepresented among outlier genes.

## Packages:
library(tidyverse)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/")

## Read genes present in the genome annotation of P. muralis (our reference genome):
ann <- read.table(file = "baypass/annotations/GCF_004329235.1_PodMur_1.0_genomic_patched_candidate_region_genes.bed", header = FALSE)
names(ann) <- c("CHROM", "START", "END", "gene", "ID")

## Restrict the universe to genes with valid gene symbols:
## This step avoids a mismatch in the universe used for the Fisher’s exact test.
## Otherwise, we might inflate the background and bias enrichment tests or distort odds ratios.
ann_sym <- ann[!is.na(ann$gene) & ann$gene != "", ]
n_total_genes <- length(unique(ann_sym$gene))
n_total_genes

## Color gene database:
ia_db <- read.csv(file = "color_genes/data/coloration_gene_database_in_Pmu_genome_annotation.csv", header = TRUE)

## For contingency table, we'll need the total number of color genes in the genome:
n_total_color_genes <- length(unique(ia_db$gene))
n_total_color_genes

## To store results of Fisher's exact test:
res_eft <- list()

## Outlier cutoff (used in BayPass and Tajima's D distributions to call outlier genes):
cutoff <- 0.9925

## Or, loop over cutoffs:
#for (cutoff in c(0.9900, 0.9925, 0.9950)) {

  ## Cutoff format:
  cutoff <- formatC(cutoff, format = "f", digits = 4)
  
  ## All outlier genes for our study species (under cutoff):
  out_genes <- read.csv(file = paste0("outliers/data_outliers/outlier_sharing_", cutoff, ".csv"), header = TRUE)
  length(out_genes$gene)
  
  ## Testing: Species:
  #sp <- "cretensis"
  
  ## Or, loop over all species:
  taxa <- sort(c("filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"))
  for (sp in taxa) {
  
    ## Outlier genes for species:
    out_sp <- out_genes$gene[out_genes[[sp]] == 1 & !is.na(out_genes$gene)]
    out_sp <- unique(out_sp)

    ## For contingency table, we'll need the total number of outlier genes for species:
    n_out_sp <- length(out_sp)
    n_out_sp
    
    ## Coloration genes among outliers of species:
    color_out_sp <- out_sp[out_sp %in% ia_db$gene]
    color_out_sp <- unique(color_out_sp)
    
    ## For contingency table, we'll need the number of outlier color genes for species:
    n_color_out_sp <- length(color_out_sp)
    n_color_out_sp
    
    ## Let's finally calculate the numbers for the contingency table:
    A <- n_color_out_sp
    B <- n_out_sp - A
    C <- n_total_color_genes - A
    D <- n_total_genes - (A + B + C)
    
    ## Compose contingency table:
    ct <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE, 
                 dimnames = list(c("outlier", "non_outlier"), c("coloration", "not_coloration")))
    
    ## Perform Fisher’s exact test (one-sided for enrichment):
    ft <- fisher.test(ct, alternative = "greater")
    
    ## Organize results:
    ft <- data.frame(cutoff = cutoff,
                     species = sp,
                     n_outlier = n_out_sp,
                     n_color_outlier = n_color_out_sp,
                     n_total_genes = n_total_genes,
                     n_total_color_genes = n_total_color_genes,
                     odds_ratio = ft$estimate,
                     p_value = ft$p.value)
    
    ## Save contingency table:
    ft$A <- A
    ft$B <- B
    ft$C <- C
    ft$D <- D

    ## Store results:
    res_eft[[length(res_eft) + 1]] <- ft

  } ## Close loop (sp).

#} ## Close loop (cutoff).

## Combine into single dataframe:
res_eft_df <- bind_rows(res_eft)

## Multiple testing correction:
res_eft_df$padj <- p.adjust(res_eft_df$p_value, method = "BH")

## Percentage of outliers corresponding to color genes:
res_eft_df$color_gene_prop <- res_eft_df$n_color_outlier / res_eft_df$n_outlier

## View:
View(res_eft_df)

## Save to file:
write.csv(res_eft_df, "color_genes/Fishers_exact_color_enrichment.csv", row.names = FALSE)

## End of script.
