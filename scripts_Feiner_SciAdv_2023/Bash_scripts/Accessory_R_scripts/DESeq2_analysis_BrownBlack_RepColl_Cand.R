dev.off()
rm(list=ls())

library(DESeq2)
library("vsn")
library("pheatmap")
library(ggplot2)
library(tidyr)
library(qqman)
library(gridExtra)
library(EnhancedVolcano)

#Read in sample sheet (metadata) and counttable
sampletable <- read.csv2("colData_all_RepCollapsed.csv", header=T, sep="," )
rownames(sampletable) <- sampletable$target_id
sampletable <- sampletable[- grep("LO102", sampletable$target_id),] #remove outlier
sampletable_BrownBlack <- sampletable[- grep("Green", sampletable$pigmentation),] #specific comparison

count_matrix_Cand <- read.csv2("raw_counts_matrix_Cand_RepCollapsed_rounded.csv", header=T, sep=",", row.names = 1)
count_matrix <- read.csv2("raw_counts_matrix_RepCollapsed_rounded.csv", header=T, sep=",", row.names = 1)
count_matrix <- count_matrix[-c(16743:16761),]
count_matrix <- rbind(count_matrix,count_matrix_Cand)
count_matrix <- count_matrix[,- grep("LO102", colnames(count_matrix))] #remove outlier
count_matrix_BrownBlack <- count_matrix[,- grep("green", colnames(count_matrix))] #specific comparison

dds_BrownBlack <- DESeqDataSetFromMatrix(countData = count_matrix_BrownBlack, colData = sampletable_BrownBlack, 
                                                                                     design = ~ pigmentation)
dds_BrownBlack <- DESeq(dds_BrownBlack)
dds_BrownBlack$pigmentation <- factor(dds_BrownBlack$pigmentation, levels=c("Brown", "Black")) #this is simply for making the plots more intuitive

### This code creates the plots in fig. S12
d <- plotCounts(dds_BrownBlack, gene="MKX", intgroup="pigmentation",returnData=TRUE)
P1 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("MKX") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="RAB18", intgroup="pigmentation",returnData=TRUE)
P2 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("RAB18") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd32", intgroup="pigmentation",returnData=TRUE)
P3 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd32") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd35", intgroup="pigmentation",returnData=TRUE)
P4 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd35") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pks1", intgroup="pigmentation",returnData=TRUE)
P5 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pks1") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pks2", intgroup="pigmentation",returnData=TRUE)
P6 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pks2") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pks3", intgroup="pigmentation",returnData=TRUE)
P7 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pks3") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd3X1", intgroup="pigmentation",returnData=TRUE)
P8 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd3X1") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd34", intgroup="pigmentation",returnData=TRUE)
P9 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd34") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pks4", intgroup="pigmentation",returnData=TRUE)
P10 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pks4") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd3", intgroup="pigmentation",returnData=TRUE)
P11 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd3") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="cyp2k6", intgroup="pigmentation",returnData=TRUE)
P12 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("cyp2k6") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ymel1l", intgroup="pigmentation",returnData=TRUE)
P13 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ymel1l") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pkssense", intgroup="pigmentation",returnData=TRUE)
P14 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pkssense") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="pks5", intgroup="pigmentation",returnData=TRUE)
P15 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("pks5") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="ptchd33", intgroup="pigmentation",returnData=TRUE)
P16 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("ptchd33") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="mastl", intgroup="pigmentation",returnData=TRUE)
P17 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("mastl") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="acbd5", intgroup="pigmentation",returnData=TRUE)
P18 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("acbd5") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="cetn2", intgroup="pigmentation",returnData=TRUE)
P19 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("cetn2") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())
d <- plotCounts(dds_BrownBlack, gene="abil", intgroup="pigmentation",returnData=TRUE)
P20 <- ggplot(d, aes(x=pigmentation, y=count, color=pigmentation)) + 
  geom_boxplot() + geom_point(position=position_jitter(), size=3) + 
  scale_color_manual(values = c("#BD912D", "black")) +
  theme_bw() + ggtitle("abil") + theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())

pdf("./Plots/DEseq_CandGenes_BrownBlack.pdf", width = 20  , height = 8, useDingbats=FALSE)
grid.arrange(P1,P2, P3, P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,
             ncol=5)
dev.off()

#get differentially expressed genes
res_BrownBlack <- results(dds_BrownBlack, contrast=c("pigmentation","Brown","Black"), alpha=0.05)
resFilt_BrownBlack <- as.data.frame(res_BrownBlack)
resFilt_BrownBlack$Genes <- rownames(resFilt_BrownBlack)
resFilt_BrownBlack_sig <- subset(resFilt_BrownBlack, padj < 0.05)
write.table(resFilt_BrownBlack_sig, "DE_Genes_BrownBlack", row.names = F, col.names = T, quote = F, sep="\t")

### This code creates Fig. 4E
#see here   https://github.com/kevinblighe/EnhancedVolcano
candidates <- c("MKX","RAB18","pks1","pks2","ptchd3","cyp2k6","ymel1l","pks3","pks4","pks5","mastl","acbd5","cetn2","abil","ptchd32","ptchd35","ptchd3X1","ptchd34","pkssense","ptchd33")
color_genes <- c('DUSP14','GNAS','PCBD1','RDH5','RDH8','RPE65')
keyvals.colour <- ifelse(rownames(res_BrownBlack) %in% candidates, "royalblue",
                         ifelse(rownames(res_BrownBlack) %in% color_genes, "black", "grey"))
#keyvals.colour <- c(rep('grey', 24927),rep('royalblue', 20))
names(keyvals.colour) <- ifelse(rownames(res_BrownBlack) %in% candidates, 'candidate_region', 
                                ifelse(rownames(res_BrownBlack) %in% color_genes,"color genes", 'genome-wide'))
keyvals.size <- ifelse(rownames(res_BrownBlack) %in% candidates, 3,
                       ifelse(rownames(res_BrownBlack) %in% color_genes,3, 2))
keyvals.alpha <- ifelse(rownames(res_BrownBlack) %in% candidates, 1,
                        ifelse(rownames(res_BrownBlack) %in% color_genes,1, 1))

P1 <- EnhancedVolcano(res_BrownBlack, lab = rownames(res_BrownBlack), x = 'log2FoldChange',y = 'padj',pCutoff = 0.05,
                      selectLab = c('DUSP14','GNAS','PCBD1','RDH5','RDH8','RPE65'), 
                      xlim=c(-8,10), ylim=c(0,17), pointSize=keyvals.size, drawConnectors = TRUE, gridlines.minor=F,gridlines.major=F,
                      widthConnectors = 0.75, lengthConnectors = unit(0.01, "npc"),shape=19,
                      colAlpha=keyvals.alpha,labSize = 6.0, colCustom = keyvals.colour)

pdf("VolcanoPlot_BrownBlack.pdf", width = 8 , height = 8, useDingbats=FALSE)
grid.arrange(P1,
             ncol=1)
dev.off()

###get positional data of genes:
Gene_positions1 <- read.table("Gene_positions_Cand", header=F)
Gene_positions2 <- read.table("Gene_positions", header=F)
Gene_positions <- rbind(Gene_positions1,Gene_positions2)
colnames(Gene_positions) <- c("Chromosome","start","end","Gene")
All_data <- merge(resFilt_BrownBlack, Gene_positions, by.x="Genes", by.y="Gene", all.x=T)
All_data$Loc <- All_data$Chromosome
All_data$Chromosome <- gsub("NC_041312.1\\b", "Chr_1", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041313.1\\b", "Chr_2", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041314.1\\b", "Chr_3", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041315.1\\b", "Chr_4", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041316.1\\b", "Chr_5", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041317.1\\b", "Chr_6", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041318.1\\b", "Chr_7", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041319.1\\b", "Chr_8", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041320.1\\b", "Chr_9", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041321.1\\b", "Chr_10", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041322.1\\b", "Chr_11", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041323.1\\b", "Chr_12", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041324.1\\b", "Chr_13", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041325.1\\b", "Chr_14", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041326.1\\b", "Chr_15", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041327.1\\b", "Chr_16", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041328.1\\b", "Chr_17", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041329.1\\b", "Chr_18", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_041330.1\\b", "Chr_19", All_data$Chromosome)
All_data$Chromosome <- gsub("NC_011607.1\\b", "Chr_20", All_data$Chromosome)
All_data$Loc <- gsub("NC_041312.1\\b", "CM014743.1", All_data$Loc)
All_data$Loc <- gsub("NC_041313.1\\b", "CM014744.1", All_data$Loc)
All_data$Loc <- gsub("NC_041314.1\\b", "CM014745.1", All_data$Loc)
All_data$Loc <- gsub("NC_041315.1\\b", "CM014746.1", All_data$Loc)
All_data$Loc <- gsub("NC_041316.1\\b", "CM014747.1", All_data$Loc)
All_data$Loc <- gsub("NC_041317.1\\b", "CM014748.1", All_data$Loc)
All_data$Loc <- gsub("NC_041318.1\\b", "CM014749.1", All_data$Loc)
All_data$Loc <- gsub("NC_041319.1\\b", "CM014750.1", All_data$Loc)
All_data$Loc <- gsub("NC_041320.1\\b", "CM014751.1", All_data$Loc)
All_data$Loc <- gsub("NC_041321.1\\b", "CM014752.1", All_data$Loc)
All_data$Loc <- gsub("NC_041322.1\\b", "CM014753.1", All_data$Loc)
All_data$Loc <- gsub("NC_041323.1\\b", "CM014754.1", All_data$Loc)
All_data$Loc <- gsub("NC_041324.1\\b", "CM014755.1", All_data$Loc)
All_data$Loc <- gsub("NC_041325.1\\b", "CM014756.1", All_data$Loc)
All_data$Loc <- gsub("NC_041326.1\\b", "CM014757.1", All_data$Loc)
All_data$Loc <- gsub("NC_041327.1\\b", "CM014758.1", All_data$Loc)
All_data$Loc <- gsub("NC_041328.1\\b", "CM014759.1", All_data$Loc)
All_data$Loc <- gsub("NC_041329.1\\b", "CM014760.1", All_data$Loc)
All_data$Loc <- gsub("NC_041330.1\\b", "CM014761.1", All_data$Loc)
All_data_sig <- subset(All_data, padj < 0.05)
write.table(All_data_sig[,c(11,9,10)], "Outliers_RNAseq_brownblack", row.names = F, col.names = F, quote = F, sep="\t")

All_data <- All_data[- grep("NW_", All_data$Chromosome),]
All_data$Log_padj <- -log10(All_data$padj)
All_data$Chromosome <- factor(All_data$Chromosome, levels = c("Chr_1","Chr_2","Chr_3","Chr_4","Chr_5","Chr_6","Chr_7","Chr_8","Chr_9","Chr_10","Chr_11","Chr_12","Chr_13","Chr_14","Chr_15","Chr_16","Chr_17","Chr_18","Chr_19","Chr_20"))
All_data <- All_data[order(All_data$start),]
All_data <- All_data[order(All_data$Chromosome),]
All_data <- separate(data = All_data, col = Chromosome, into = c("scaffold_name", "scaffold_Nr"), sep = "_")
All_data$scaffold_Nr <- as.numeric(All_data$scaffold_Nr)
pdf("DEseq_manhattan_BrownBlack.pdf", width = 14  , height = 8, useDingbats=FALSE)
manhattan(All_data, snp="Genes", chr="scaffold_Nr", bp="start",  p="padj", logp = T, suggestiveline = F, genomewideline = F, col = c("#FF0073", "#8cece2"))
dev.off()
pdf("DEseq_manhattan_BrownBlack_Chr12.pdf", width = 14  , height = 8, useDingbats=FALSE)
manhattan(subset(All_data,scaffold_Nr==12), snp="Genes", chr="scaffold_Nr", bp="start",  p="padj", logp = T, suggestiveline = F, genomewideline = F, col = c("#FF0073", "#8cece2"))
dev.off()
All_data_Chr12 <- subset(All_data, scaffold_Nr == "12")
All_data_Chr12$Log_p <- -log10(All_data_Chr12$pvalue)
All_data_Chr12$Sign <- "No"
All_data_Chr12$Sign[All_data_Chr12$padj < 0.01] <- "Yes"
p1 <- ggplot(All_data_Chr12, aes(x=start, y=Log_p)) +
  geom_point(aes(color=Sign),size=2) + scale_color_manual(values = c("brown","red")) +
  xlim(21300000, 21900000) +
  theme_bw() + theme(legend.position="none", plot.title = element_text(size=20)) 
pdf("DEseq_manhattan_BrownBlack_cand_RepColl.pdf", width = 14  , height = 8, useDingbats=FALSE)
grid.arrange(p1, ncol=1)
dev.off()

write.csv(All_data, file="resFilt_BrownBlack_pos.csv")
