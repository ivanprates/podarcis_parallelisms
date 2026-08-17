dev.off()
rm(list=ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)
library(dplyr)
options(scipen=999)

#do this for both IT and SA
Data_brown <- read.table("WG_5kb_TajimaD_brown_IT.Tajima.D", header = T)
Data_green <- read.table("WG_5kb_TajimaD_green_IT.Tajima.D", header = T)
Data_brown$POS <- paste(Data_brown$CHROM,Data_brown$BIN_START, sep="_")
Data_green$POS <- paste(Data_green$CHROM,Data_green$BIN_START, sep="_")
Data <- merge(Data_green, Data_brown, by.x="POS", by.y="POS")
Data <- Data[,-c(2,3,4,6:8)]
colnames(Data) <- c("POS", "green", "brown")
Data$DeltaTajima <- Data$brown-Data$green
Data[c('CHROM', 'BP')] <- str_split_fixed(Data$POS, '_', 2)
Data <- Data[- grep("NaN", Data$DeltaTajima),]
upp_quant <- quantile(Data$DeltaTajima, 0.9975, na.rm = T)
Data <- Data %>% mutate(outlier = ifelse(DeltaTajima > upp_quant, "outlier", "background"))
Data %>% group_by(outlier) %>% tally()
Data_outlier <- subset(Data, outlier == "outlier")
Data_outlier <- Data_outlier[,c(5,6)]
Data_outlier$end <- as.numeric(Data_outlier$BP)+5000
colnames(Data_outlier) <- c("Loc", "start","end")
Data_outlier <- Data_outlier[order(as.numeric(Data_outlier$start)),]
Data_outlier <- Data_outlier[order(Data_outlier$Loc),]
write.table(Data_outlier,"Outliers_TajimasD_IT.bed", sep = "\t", col.names=F, row.names=F, quote=F)
