### This code generates fig. S12a-c

dev.off()
rm(list=ls())

library(DESeq2)
library("vsn")
library("pheatmap")
library(ggplot2)
library(tidyr)
library("RColorBrewer")

#Read in sample sheet (metadata) and counttable
sampletable <- read.csv2("colData_all.csv", header=T, sep=",")
rownames(sampletable) <- sampletable$target_id

count_matrix_Cand <- read.delim("raw_counts_matrix_Cand.txt", header=T, sep="\t", row.names = 1)
count_matrix <- read.delim("raw_counts_matrix.txt", header=T, sep="\t", row.names = 1)
count_matrix <- count_matrix[-c(16743:16761),]
count_matrix <- rbind(count_matrix,count_matrix_Cand)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sampletable, design = ~ lineage + pigmentation)
dds <- DESeq(dds)

###Variance stabalization for Heatmap, PCA, etc.
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("pigmentation", "individual"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./Plots/DEseq_PCA_AllSamples.pdf", width = 5  , height = 5, useDingbats=FALSE)
ggplot(pcaData, aes(PC1, PC2, shape=individual, color=pigmentation)) +  theme_bw() +
  geom_point(size=3) + scale_shape_manual(values=1:10) +
  scale_colour_manual(values = c("black", "#BD912D", "#72A042")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#Same with LO102 excluded
sampletable <- sampletable[- grep("LO102", sampletable$target_id),] #remove outlier
count_matrix <- count_matrix[,- grep("LO102", colnames(count_matrix))] #remove outlier
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sampletable, design = ~ lineage + pigmentation)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("pigmentation", "individual","lineage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./Plots/DEseq_PCA_LO102excl_col.pdf", width = 5  , height = 5, useDingbats=FALSE)
ggplot(pcaData, aes(PC1, PC2, shape=individual, color=pigmentation)) +  theme_bw() +
  geom_point(size=3) + scale_shape_manual(values=c(1:5,7:10)) +
  scale_colour_manual(values = c("black", "#BD912D", "#72A042")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio=1) +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#Make heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$pigmentation, vsd$ind.n, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(c("white","lightgrey","grey80","turquoise","deepskyblue3","blue3","navyblue","black")))(255)
pdf("./Plots/DEseq_heatmap_LO102excl.pdf", width = 10  , height = 10, useDingbats=FALSE)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
