rm(list=ls())
dev.off()

library("genoPlotR")
library(gridExtra)


setwd("C:/Users/Nathalie/Dropbox/WallLizard_NanoporeAssemblies/")

############### Fig. 2D - gene annotations in candidate region (exon level)

data <- read.csv("CandidateRegion_Genes_Location_Spa_exons.csv", header = T)
dna_data <- dna_seg(data)
dna_seg <- list(dna_data)
annot1 <- annotation(x1=middle(dna_data), text=dna_data$name, rot=65)
genes_Pod <- unique(data$name)
x_Pod <- sapply(genes_Pod, function(x)
  range(dna_data
        [dna_data$name == x,]))
annot_Pod <- annotation(x1=x_Pod[1,], x2=x_Pod[2,], text=dimnames(x_Pod)[[2]], rot=30) 
pdf("./Cand_region_genes_exons.pdf", width = 14  , height = 21.5, useDingbats=FALSE)
plot_gene_map(dna_segs=dna_seg, gene_type = "side_blocks", annotations = annot_Pod,annotation_height=3,
              dna_seg_scale=TRUE, seg_plot_height=10)
dev.off()

############### Fig. 3B - gene annotations in candidate region (gene level), comparison between different genomes
                
pdf("./Cand_region_genes_comp.pdf", width = 14  , height = 21.5, useDingbats=FALSE)
data <- read.csv("CandidateRegion_Genes_Location_Spa.csv", header = T)
dna_data <- dna_seg(data)
dna_seg <- list(dna_data)
annot1 <- annotation(x1=middle(dna_data), text=dna_data$name, rot=65)
plot_gene_map(dna_segs=dna_seg, scale = F, gene_type = "side_blocks", annotations = annot1, seg_plot_height=10, 
              dna_seg_scale=TRUE, main="Spanish reference genome")

data <- read.csv("CandidateRegion_Genes_Location_SMA_HiC_reverse.csv", header = T)
dna_data <- dna_seg(data)
dna_seg <- list(dna_data)
annot1 <- annotation(x1=middle(dna_data), text=dna_data$name, rot=65)
xlims <- list(c(40999000,40560000)) 
plot_gene_map(dna_segs=dna_seg, scale = F, gene_type = "side_blocks", annotations = annot1, seg_plot_height=10, 
                    dna_seg_scale=TRUE, xlims = xlims, main="SMA_HiC") 

data <- read.csv("CandidateRegion_Genes_Location_V1_CDS_HiC.csv", header = T)
dna_data <- dna_seg(data)
dna_seg <- list(dna_data)
annot1 <- annotation(x1=middle(dna_data), text=dna_data$name, rot=65)
plot_gene_map(dna_segs=dna_seg, scale = F, gene_type = "side_blocks", annotations = annot1, seg_plot_height=10, 
                    dna_seg_scale=TRUE, main="CDS_HiC") 

dev.off()
