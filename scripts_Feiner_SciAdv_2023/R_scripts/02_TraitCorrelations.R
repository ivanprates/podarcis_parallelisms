rm(list=ls())

library(openxlsx)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(pcaMethods)
library(dplyr)

options(scipen=999)

#Import the current new data set - this includes data up till 2018
data <- read.csv("Phenotype_data_final_June2023.csv")

#Import lineage assignments per population (from structure on RAD-seq data plus extrapolation)
Pops <- read.csv("PopsLocs.csv")
Pops <- Pops[,-c(3,4)]
data <- merge(data, Pops, by = "abbpop", all.x=T)

#Greenness per pop
Green_mean_sex <- as.data.frame(data %>% group_by(abbpop, sex) %>% dplyr::summarize(Green = mean(GreenResolved, na.rm=TRUE)))
Green_mean_sex <- subset(Green_mean_sex, sex=="M")
Green_mean_sex <- Green_mean_sex[,-c(2)]
#Green_mean <- data %>% group_by(abbpop) %>% dplyr::summarize(Green = mean(GreenResolved, na.rm=TRUE))
data <- merge(data, Green_mean_sex, by="abbpop")

Pops <- merge(Pops, Green_mean_sex, by="abbpop")
Pops_SA_intro <- subset(Pops, lineage=="SA" & Green>=3)
Pops_SA_intro$lineage_intro <- "SA_intro"
Pops_SA_intro <- Pops_SA_intro[,-c(2,3)]

data <- merge(data,Pops_SA_intro, by="abbpop", all.x=T)

#change some colnames
colnames(data)[c(11,14,15,16)] <- c("green","black","blue","rel.headlength")

#subset data by lineage
data_IT_M <- subset(data, lineage=="ITA" & sex=="M")
data_Hyb_M <- subset(data, lineage=="MIXED" & sex=="M")
data_SAintro_M <- subset(data, lineage_intro == "SA_intro" & sex=="M")
data_IT_F <- subset(data, lineage=="ITA" & sex=="F")
data_Hyb_F <- subset(data, lineage=="MIXED" & sex=="F")
data_SAintro_F <- subset(data, lineage_intro == "SA_intro" & sex=="F")

############### Fig. 1D - pairwise correlation heatmaps

#correlation IT males
matrix_IT_M <- as.matrix(scale(data_IT_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("IT_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation hybrid males
matrix_IT_M <- as.matrix(scale(data_Hyb_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("Hyb_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation SAintro males
matrix_IT_M <- as.matrix(scale(data_SAintro_M[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_M, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("SAintro_M_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation IT females
matrix_IT_F <- as.matrix(scale(data_IT_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("IT_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation hybrid females
matrix_Mix_F <- as.matrix(scale(data_Hyb_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_Mix_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("Hyb_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

#correlation SAintro females
matrix_IT_F <- as.matrix(scale(data_SAintro_F[,c(8,16,11,14,15)]))
cormat <- rcorr(matrix_IT_F, type = "pearson")
cormat_cor <- cormat$r
cormat_p <- cormat$P
pdf("SAintro_F_corrplot.pdf", height = 5,  useDingbats=FALSE)
corrplot(cormat_cor, p.mat = cormat_p, tl.col ="black", method = 'color',  type = 'upper', order="original",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.4, insig = 'label_sig', pch.col = 'grey20', tl.srt = 45)
dev.off()

###Phenotypic integration score
library(geomorph)

#compare integration in males between IT, Hyb and SAintro
data_M <- subset(data, (lineage_intro=="SA_intro" | lineage== "ITA" | lineage== "MIXED") & sex=="M")

#delete incomplete obs
data_M <- data_M[complete.cases(data_M[,c(8,16,11,14,15)]), ]

#scale the data and split by lineage
data_M[,c(8,16,11,14,15)] <- scale(data_M[,c(8,16,11,14,15)])
data_M_IT <- subset(data_M, lineage=="ITA")
data_M_Hyb <- subset(data_M, lineage=="MIXED")
data_M_SAintro <- subset(data_M, lineage_intro=="SA_intro")

IT_M <- integration.Vrel(as.matrix(data_M_IT[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
HY_M <- integration.Vrel(as.matrix(data_M_Hyb[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
SA_M <- integration.Vrel(as.matrix(data_M_SAintro[,c(8,16,11,14,15)])) #extract $Re.obs as Vrel
res <- compare.ZVrel(IT_M, HY_M, SA_M) #extract P-value for difference

#same for females
data_F <- subset(data, (lineage_intro=="SA_intro" | lineage== "ITA" | lineage== "MIXED") & sex=="F")

data_F <- data_F[complete.cases(data_F[,c(8,16,11,14,15)]), ]

data_F[,c(8,16,11,14,15)] <- scale(data_F[,c(8,16,11,14,15)])
data_F_IT <- subset(data_F, lineage=="ITA")
data_F_Hyb <- subset(data_F, lineage=="MIXED")
data_F_SAintro <- subset(data_F, lineage_intro=="SA_intro")

IT_F <- integration.Vrel(as.matrix(data_F_IT[,c(8,16,11,14,15)]))
HY_F <- integration.Vrel(as.matrix(data_F_Hyb[,c(8,16,11,14,15)]))
SA_F <- integration.Vrel(as.matrix(data_F_SAintro[,c(8,16,11,14,15)]))
res <- compare.ZVrel(IT_F, HY_F, SA_F)
