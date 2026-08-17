### This code creates fig. S13

dev.off()
rm(list=ls())

library(tidyverse)
library(gridExtra)
library(dplyr)
library(Hmisc)

options(scipen=999)

Data_chr1 <- read.table("LD_decay_chr1.ld_decay", header=T)
Data_chr2 <- read.table("LD_decay_chr2.ld_decay", header=T)
Data_chr3 <- read.table("LD_decay_chr3.ld_decay", header=T)
Data_chr4 <- read.table("LD_decay_chr4.ld_decay", header=T)
Data_chr5 <- read.table("LD_decay_chr5.ld_decay", header=T)
Data_chr6 <- read.table("LD_decay_chr6.ld_decay", header=T)
Data_chr7 <- read.table("LD_decay_chr7.ld_decay", header=T)
Data_chr8 <- read.table("LD_decay_chr8.ld_decay", header=T)
Data_chr9 <- read.table("LD_decay_chr9.ld_decay", header=T)
Data_chr10 <- read.table("LD_decay_chr10.ld_decay", header=T)
Data_chr11 <- read.table("LD_decay_chr11.ld_decay", header=T)
Data_chr12 <- read.table("LD_decay_chr12.ld_decay", header=T)
Data_chr13 <- read.table("LD_decay_chr13.ld_decay", header=T)
Data_chr14 <- read.table("LD_decay_chr14.ld_decay", header=T)
Data_chr15 <- read.table("LD_decay_chr15.ld_decay", header=T)
Data_chr16 <- read.table("LD_decay_chr16.ld_decay", header=T)
Data_chr17 <- read.table("LD_decay_chr17.ld_decay", header=T)
Data_chr18 <- read.table("LD_decay_chr18.ld_decay", header=T)
Data_chrZ <- read.table("LD_decay_chrZ.ld_decay", header=T)

Data <- rbind(Data_chr1,Data_chr2,Data_chr3,Data_chr4,Data_chr5,Data_chr6,Data_chr7,Data_chr8,Data_chr9,Data_chr10,
              Data_chr11,Data_chr12,Data_chr13,Data_chr14,Data_chr15,Data_chr16,Data_chr17,Data_chr18,Data_chrZ)

Data$chr <- gsub("CM014743.1\\b", "Chr 1", Data$chr)
Data$chr <- gsub("CM014744.1\\b", "Chr 2", Data$chr)
Data$chr <- gsub("CM014745.1\\b", "Chr 3", Data$chr)
Data$chr <- gsub("CM014746.1\\b", "Chr 4", Data$chr)
Data$chr <- gsub("CM014747.1\\b", "Chr 5", Data$chr)
Data$chr <- gsub("CM014748.1\\b", "Chr 6", Data$chr)
Data$chr <- gsub("CM014749.1\\b", "Chr 7", Data$chr)
Data$chr <- gsub("CM014750.1\\b", "Chr 8", Data$chr)
Data$chr <- gsub("CM014751.1\\b", "Chr 9", Data$chr)
Data$chr <- gsub("CM014752.1\\b", "Chr 10", Data$chr)
Data$chr <- gsub("CM014753.1\\b", "Chr 11", Data$chr)
Data$chr <- gsub("CM014754.1\\b", "Chr 12", Data$chr)
Data$chr <- gsub("CM014755.1\\b", "Chr 13", Data$chr)
Data$chr <- gsub("CM014756.1\\b", "Chr 14", Data$chr)
Data$chr <- gsub("CM014757.1\\b", "Chr 15", Data$chr)
Data$chr <- gsub("CM014758.1\\b", "Chr 16", Data$chr)
Data$chr <- gsub("CM014759.1\\b", "Chr 17", Data$chr)
Data$chr <- gsub("CM014760.1\\b", "Chr 18", Data$chr)
Data$chr <- gsub("CM014761.1\\b", "Chr Z", Data$chr)
Data$chr <- factor(Data$chr, levels = c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7","Chr 8","Chr 9","Chr 10","Chr 11","Chr 12","Chr 13","Chr 14","Chr 15","Chr 16","Chr 17","Chr 18","Chr Z"))

pdf("LD_decay_all_chromosomes.pdf", width = 9, height = 7, useDingbats=FALSE)
jpeg("LD_decay_all_chromosomes.jpg", width = 1800, height = 1400, res=300)

ggplot(Data, aes(distance, avg_R2)) + xlim(1,10000)+ ylim(0,0.4) +
  geom_jitter(alpha=0.2, size=0.1) +
  geom_smooth(se=F, aes(color=chr)) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) + theme_bw()
dev.off()
