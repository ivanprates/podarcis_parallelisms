## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2025.
## The goals of this R script are:
## To plot scans of genome-wide differentiation from BayPass results.
## To plot patterns of outlier gene sharing across Podarcis lizard species.
## To compose manuscript figures.

## 1. Getting ready ----

## Packages:
library(tidyverse)

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Path to working directory: Ivan's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## 2. Outlier genes ----

## List all files with unique outlier genes:
files <- list.files(path = "outliers/data_outliers", 
                    pattern = "_SNPs_in_TajimasD_windows_10kb_outlier_genes_unique_.*\\.txt$", 
                    full.names = TRUE)

## Read and combine:
g_ab <- do.call(rbind, lapply(files, function(f) {
  dat <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(dat) <- c("CHROM", "gene", "ID")
  dat$source_file <- basename(f)
  dat$source_file <- gsub(x = dat$source_file, pattern = "_SNPs_in_TajimasD_windows_10kb_outlier_genes_unique|\\.txt", replacement = "")
  dat <- dat %>% separate(source_file, into = c("species", "cutoff"), sep = "_", remove = TRUE)
}))

## Capitalize residual gene symbols:
g_ab$gene <- str_to_upper(g_ab$gene)

## Number of outlier genes for cutoff combination:
gene_counts <- g_ab %>%
  group_by(species, cutoff) %>%
  summarise(N_IDs = n_distinct(ID), .groups = "drop")

## As numeric:
gene_counts$cutoff <-as.numeric(gene_counts$cutoff)

## Plot:
pg <- ggplot(data = gene_counts, aes(x = cutoff, y = N_IDs, color = species, group = species)) +
  geom_line() +
  geom_point() +
  labs(x = "Cutoff (for both Tajima's D and BayPass SNPs distributions)", 
       y = "Number of unique candidate genes",
       color = "Species") +
  scale_x_continuous(breaks = sort(unique(gene_counts$cutoff))) +
  theme_bw()

## Save::
ggsave(plot = pg, width = 20, height = 15, units = "cm", dpi = 300, 
       filename = "TajimasD/0326/plots/selection_of_cutoffs.jpg")

## Select cutoff:
#c <- 0.9925

## Or, loop over combinations:
for (c in c(0.9900, 0.9925, 0.9950)) {

## Format:
c <- formatC(c, format = "f", digits = 4)

## Subset data:
g_df <- g_ab[g_ab$cutoff == c, c("ID", "gene", "species")]
g_df$present <- 1

## Remove residual duplicated gene symbols (from manual patching of the candidate region):
g_df <- g_df[!duplicated(g_df[, c("gene", "species")]), ]

## To wide format:
g_df <- g_df %>% 
  pivot_wider(
    id_cols = c(ID, gene),               
    names_from = species,                
    values_from = present,               
    values_fill = list(present = 0))

## Column names:
taxa <- sort(c("filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae"))
names(g_df) <- c("ID", "gene", taxa)

## Sort:
g_df <- arrange(g_df, by = gene)

## Sum across columns for each row to get number of outlier windows shared across taxa:
g_df$shared_by <- rowSums(g_df[taxa], na.rm = TRUE)

## Import P. muralis candidate genes from Feiner et al. 2024 Sci Adv ----
Pmu_genes <- read.table(file = "baypass/annotations/Feiner_2024_Outliers_IT_SA_5kb2.5_FST_merged_d10kb_10kb_up_and_down_outlier_genes_unique.txt", header = FALSE)
names(Pmu_genes) <- c("CHROM", "gene", "ID")
Pmu_genes$gene <- str_to_upper(Pmu_genes$gene)
Pmu_genes <- Pmu_genes[!duplicated(Pmu_genes$gene), ] ## Remove lingering duplicates, if they exist. 
Pmu_genes$Pmu_candidate <- 1

## Merge:
g_sh <- merge(g_df, Pmu_genes[c("ID", "Pmu_candidate")], by = "ID", all.x = TRUE)
g_sh$Pmu_candidate[is.na(g_sh$Pmu_candidate)] <- 0

## Arrange:
g_sh <- arrange(g_sh, ID)
g_sh <- arrange(g_sh, desc(shared_by))

## Save:
write.csv(g_sh, 
          file = paste0("TajimasD/0326/data_outliers/outlier_sharing_", c, ".csv"), quote = FALSE, row.names = FALSE)

## Some numbers:
## How many genes in species sets with different numbers of species?
table(rowSums(g_sh[taxa], na.rm = TRUE))

## How many outlier genes per species?
apply(X = g_sh[taxa], MARGIN = 2, FUN = sum)

} ## Close loop (c).

## End of script.
