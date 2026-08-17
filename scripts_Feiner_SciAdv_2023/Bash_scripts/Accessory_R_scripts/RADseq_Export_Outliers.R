dev.off()
rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(gridExtra)
library(dplyr)
library(qqman)
options(scipen=999)

Plink_adjusted <- read.table("plink.assoc.linear.adjusted", header=T) #do this for both IT and SA
Plink_nonadjusted <- read.table("plink.assoc.linear", header=T) #do this for both IT and SA
Plink_nonadjusted <- Plink_nonadjusted[,1:3]
Plink_nonadjusted <- Plink_nonadjusted[!duplicated(Plink_nonadjusted), ]
Plink_new <- merge(Plink_adjusted, Plink_nonadjusted, by.x = "SNP", by.y = "SNP", all.x=T, sort=F)
Plink_new <- Plink_new[,c(2,1,12,3:10)]
quantile(Plink_new$FDR_BH, c(0.975, 0.995), na.rm = T)
Plink_new <- Plink_new %>% mutate(outlier = ifelse(FDR_BH < 0.05, "outlier", "background"))
Plink_new %>% group_by(outlier) %>% tally()
Data_outlier <- subset(Plink_new, outlier == "outlier")
Data_outlier <- Data_outlier[,-c(2,4:9,11,12)]
Data_outlier <- Data_outlier[order(Data_outlier$BP),]
Data_outlier <- Data_outlier[order(Data_outlier$CHR.x),]
Data_outlier$start <- Data_outlier$BP-5000 #test with different window sizes (5000 and 50000)
Data_outlier$end <- Data_outlier$BP+5000 #test with different window sizes (5000 and 50000)
Data_outlier$start[Data_outlier$start<0] <- 1
Data_outlier <- Data_outlier[,c(1,4,5)]
write.table(Data_outlier,"Outliers_RAD_IT_new.bed", sep = "\t", col.names=F, row.names=F, quote=F)
