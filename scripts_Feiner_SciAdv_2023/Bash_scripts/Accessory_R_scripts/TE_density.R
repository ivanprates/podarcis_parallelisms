### This code creates fig. S6a

dev.off()
rm(list=ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)

Data <- read.table("Chr12_binning_summary_10kb", header = F)
Data <- Data[-6111,]
Data$V2 <- as.integer(Data$V2)

pdf("TE_density.pdf", width=5, height=5, useDingbats=FALSE)
ggplot(Data, aes(x=V2, y=V1)) + geom_point(stat = 'identity') + xlab("Position on Chr 12") + ylab("Number of TEs per 10kb") +
  xlim(21000000,22000000) + ylim(0,15) + theme_bw()
dev.off()
