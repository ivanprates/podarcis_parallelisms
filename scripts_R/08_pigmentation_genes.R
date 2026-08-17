## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To compile pigment- or color-related genes from the literature, focusing on chromatophore-bearing vertebrates.
## To extract such genes that were inferred as outliers between populations of each Podarcis species.

## 1. Getting ready and some checks ----

## Packages:
library(tidyverse)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/")

## Create directories to save outputs:
dir.create("color_genes/data")

## Read genes present in the genome annotation of P. muralis:
ann <- read.table(file = "baypass/annotations/GCF_004329235.1_PodMur_1.0_genomic_patched_candidate_region_genes.bed", header = FALSE)
names(ann) <- c("CHROM", "START", "END", "gene", "ID")

## Read candidate P. muralis genes in Feiner et al. 2024:
P <- read.csv(file = "baypass/annotations/FST_outliers_overlap_IT_SA_89genes.csv", header = TRUE)
Pmu <- sort(P$external_gene_name)

## Total candidate genes in P. muralis:
length(Pmu)

## Total genes with gene symbol annotations:
Pmu <- setdiff(Pmu, "")
length(Pmu)

## Which candidate P. muralis genes having symbols are present in the annotation used?
intersect(Pmu, ann$gene)

## Not in annotation:
nia <- Pmu[!Pmu %in% ann$gene]
nia

## Read candidate region genes in the chromosome 12 of P. muralis, from Feiner et al. 2024:
R <- read.csv(file = "baypass/annotations/CandidateRegion_Genes_Location_Spa_exons.csv", header = TRUE)
R <- sort(unique(R$name))
R <- toupper(R)

## Are candidate region genes (chrom 12) present in the annotation used?
length(R)
intersect(R, ann$gene)

## Not in annotation:
R[!R %in% ann$gene]
ann[grepl(x = ann$gene, pattern = "PKS"), ]
ann[grepl(x = ann$gene, pattern = "PTCHD"), ]

## 2. Coloration gene databases ----

## Read and edit database: McLean et al. 2017:
## Revealing the biochemical and genetic basis of color variation in a polymorphic lizard.
## Mol. Biol. Evol. (2017), 34(8): 1924–1935
M <- read.csv(file = "color_genes/databases/McLean2017_MBE_ColorGenes.csv", header = TRUE)
M <- M[, -2]
names(M) <- c("gene", "pheno_McLean2017")
M$in_McLean2017 <- 1

## Column names:
names(M)[names(M) == "pheno_McLean2017"] <- "phenotype_notes_McLean2017"

## Are M genes present in the annotation used?
length(M$gene)
length(intersect(M$gene, ann$gene))
setdiff(M$gene, M$gene[M$gene %in% ann$gene]) ## Not in annotation.

## Read and edit database: Baxter et al. 2018:
## A curated gene list for expanding the horizons of pigmentation biology.
## Pigment Cell Melanoma Res. (2019), 32: 348–358.
B <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0007-tables7.csv", header = TRUE)
names(B)[names(B) == "Human.gene.symbol"] <- "gene"
names(B)[names(B) == "Gene.stable.ID"] <- "ID"
B$gene[B$gene == ""] <- B$Mouse.gene.symbol[B$gene == ""]
B$gene[B$gene == ""] <- B$Zebrafish.gene.symbol1[B$gene == ""]
B$gene <- toupper(B$gene)

## Adding phenotypes:
ph <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0001-tables1.csv", header = TRUE)
B$pheno_Baxter2018_human[B$gene %in% ph$Gene.symbol] <- ph$Phenotype.description.notes[B$gene %in% ph$Gene.symbol]

pm <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0002-tables2.csv", header = TRUE)
B$pheno_Baxter2018_mouse[B$Mouse.gene.symbol %in% pm$Mouse.gene.symbol] <- pm$Phenotype.description.notes[B$Mouse.gene.symbol %in% pm$Mouse.gene.symbol]

pz <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0003-tables3.csv", header = TRUE)
B$pheno_Baxter2018_zebrafish[B$Zebrafish.gene.symbol1 %in% pz$Zebrafish.gene.symbol] <- pz$Phenotype.description.notes[B$Zebrafish.gene.symbol1 %in% pz$Zebrafish.gene.symbol]

pg <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0004-tables4.csv", header = TRUE)
B$pheno_Baxter2018_GO[B$ID %in% pg$Gene.stable.ID] <- pg$Phenotype.descriptions.notes..for.genes.retrieved.exclusively.from.GO...phenotype.descriptions.for.all.other.genes.are.located.in.tables.for.OMIM..MGI..and.ZFIN.[B$ID %in% pg$Gene.stable.ID]

pp <- read.csv(file = "color_genes/databases/Baxter2018_pcmr12743-sup-0005-tables5.csv", header = TRUE)
B$pheno_Baxter2018_PubMed[B$ID %in% pp$Gene.stable.ID] <- pp$Phenotype.description.notes[B$ID %in% pp$Gene.stable.ID]

## Columns with slashes:
B$gene[grepl(x = B$gene, pattern = "/")]

## Split rows where IDs contain slashes:
B <- B %>% separate_rows(gene, sep = " / ")

## Exclude data on humans and mouse, keep zebrafish:
B <- B[grepl(x = B$Species.with.phenotype, pattern = "Zebrafish|zebrafish"), ]

## Simplify, annotate:
B <- B[c("gene", "pheno_Baxter2018_zebrafish")]
## Removed: "pheno_Baxter2018_human", "pheno_Baxter2018_mouse", "pheno_Baxter2018_GO", "pheno_Baxter2018_PubMed"
B$in_Baxter2018 <- 1

## Column names:
names(B)[names(B) == "pheno_Baxter2018_zebrafish"] <- "phenotype_notes_Baxter2018"

## Are B genes present in the annotation used?
length(B$gene)
length(intersect(B$gene, ann$gene))
setdiff(B$gene, B$gene[B$gene %in% ann$gene]) ## Not in annotation.

## Read and edit database: Elkin et al. 2023 (GepheBase):
## Analysis of the genetic loci of pigment pattern evolution in vertebrates.
## Biol. Rev. (2023), 98: 1250–1277.
E <- read.csv(file = "color_genes/databases/Gephebase-entries-2025-09-22.csv", header = TRUE)
E$Generic.Gene.Name <- toupper(E$Generic.Gene.Name)

## Exclude taxa:
unique(sort(E$Taxon.A.Lineage))
E <- E[!grepl(x = E$Taxon.A.Lineage, pattern = "Arthropoda"), ]
E <- E[!grepl(x = E$Taxon.A.Lineage, pattern = "Viridiplantae"), ]
E <- E[!grepl(x = E$Taxon.A.Lineage, pattern = "Aves"), ]
E <- E[!grepl(x = E$Taxon.A.Lineage, pattern = "Mammalia"), ]

## Simplify, annotate:            
E <- E[c("Generic.Gene.Name", "Trait")]
names(E) <- c("gene", "phenotype_notes_Elkin2023")
E$in_Elkin2023 <- 1
E <- E[!duplicated(E$gene), ]
E <- E[!E$gene == "", ]

## Are E genes present in the annotation used?
length(E$gene)
length(intersect(E$gene, ann$gene))
setdiff(E$gene, E$gene[E$gene %in% ann$gene]) ## Not in annotation.

## Combine databases:
color_db <- merge(M, B, by = "gene", all = TRUE)
color_db <- merge(color_db, E, by = "gene", all = TRUE)

## NAs to zero:
color_db$in_Baxter2018[is.na(color_db$in_Baxter2018)] <- 0
color_db$in_Elkin2023[is.na(color_db$in_Elkin2023)] <- 0
color_db$in_McLean2017[is.na(color_db$in_McLean2017)] <- 0

## Column order:
color_db <- color_db[c(1, 5, 7, 3, 4, 6, 2)]
dim(color_db)

## Save:
write.csv(color_db, file = "color_genes/data/coloration_gene_database.csv", row.names = FALSE)

## Combined color genes present in Pmu genome annotation:
length(color_db$gene)
length(intersect(color_db$gene, ann$gene))
ia_db <- color_db[color_db$gene %in% ann$gene, ]

## Save:
write.csv(ia_db, file = "color_genes/data/coloration_gene_database_in_Pmu_genome_annotation.csv", row.names = FALSE)

## Genes not in Pmu annotation:
nia <- setdiff(color_db$gene, color_db$gene[color_db$gene %in% ann$gene])
nia_db <- color_db[color_db$gene %in% nia, ]

## Save:
write.csv(nia_db, file = "color_genes/data/coloration_gene_database_not_in_Pmu_genome_annotation.csv", row.names = FALSE)

## 3. Coloration genes in shared outlier gene set ----

## Cutoff:
cutoff <- 0.9925

## Or, loop over combinations:
#for (c in c(0.9900, 0.9925, 0.9950)) {

  ## Format:
  cutoff <- formatC(cutoff, format = "f", digits = 4)
  
  ## Read shared outlier data:
  genes <- read.csv(file = paste0("outliers/data_outliers/outlier_sharing_", cutoff, ".csv"), header = TRUE)
  names(genes)[names(genes) == "Pmu_candidate"] <- "muralis_candidate"
  
  ## Combine with color database:
  m <- merge(genes, color_db, by = "gene", all.x = TRUE)
  
  ## NAs to zero after merging:
  m$in_Baxter2018[is.na(m$in_Baxter2018)] <- 0
  m$in_Elkin2023[is.na(m$in_Elkin2023)] <- 0
  m$in_McLean2017[is.na(m$in_McLean2017)] <- 0

  ## How many in each database:
  sum(m$in_McLean2017)
  sum(m$in_Baxter2018)
  sum(m$in_Elkin2023)
  
  ## Present in at least one pigmentation gene database:
  C <- m[m$in_McLean2017 != 0 | m$in_Baxter2018 != 0 | m$in_Elkin2023 != 0, ]
  
  ## Arrange:
  C <- arrange(C, desc(shared_by))
  
  ## Save:
  write.csv(C, file = paste0("color_genes/data/coloration_genes_", cutoff, ".csv"), row.names = FALSE)

#} ## Close loop (c).

## End of script.
