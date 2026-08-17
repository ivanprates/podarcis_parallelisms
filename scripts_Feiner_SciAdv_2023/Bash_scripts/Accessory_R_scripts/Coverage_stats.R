### This code creates Fig. 3A and fig. S5

dev.off()
rm(list=ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)

Data_depth_ind <- read.table("Read_Depth_Per_Ind", header = T, sep="\t") #can also be done with Read_Depth_Per_Ind_post
Data_depth_ind[,c(2:150)] <- sapply(Data_depth_ind[,c(2:150)], as.integer)

Data_depth_green <- Data_depth_ind[,c(1,26:32,117:119,33:35,124:126,130:133,14:17,68:73)]
Data_depth_brown <- Data_depth_ind[,c(1:3,18:25,127:129,65:67,120:123,4:13)]
Data_depth_green$sum <- rowSums(Data_depth_green[, c(2:31)], na.rm=T)
Data_depth_green <- Data_depth_green[,c(1,32)]
Data_depth_brown$sum <- rowSums(Data_depth_brown[, c(2:31)], na.rm=T)
Data_depth_brown <- Data_depth_brown[,c(1,32)]
Data_depth_ind$sum <- rowSums(Data_depth_ind[, c(2:150)], na.rm=T)
Data_depth_ind <- Data_depth_ind[,c(1,151)]
Data <- rbind(Data_depth_green,Data_depth_brown)
Data$group <- c(rep("green",28270),rep("brown",28270))
Data$group <- c(rep("green",4821),rep("brown",4821))

ggplot(Data, aes(x=POS, y=sum, col=group)) + geom_line(stat = 'identity', alpha=1) +
  scale_colour_manual(values=c("#816830","#45A848")) + ylim(0,1500) + theme_bw()

p0 <- ggplot(Data_depth_ind, aes(x=POS, y=sum)) + geom_line(stat = 'identity', col="black") + theme(aspect.ratio = 0.1) +
  theme_bw()
p1 <- ggplot(Data_depth_ind, aes(x=POS, y=sum)) + geom_line(stat = 'identity', col="black") + theme(aspect.ratio = 0.1) +
  ylim(0,6000) + theme_bw()
p2 <- ggplot(Data_depth_green, aes(x=POS, y=sum)) + geom_line(stat = 'identity', col="#45A848") + theme(aspect.ratio = 0.1) +
  ylim(0,1210) + theme_bw()
p3 <- ggplot(Data_depth_brown, aes(x=POS, y=sum)) + geom_line(stat = 'identity', col="#614E23") + theme(aspect.ratio = 0.1) +
  ylim(0,1210) + theme_bw()

#Number of SNPs per 10kb-window
Data_pre_20mb <- read.table("GATK_all_SNP_filtered_cand_plusminus20mb_coverage_10kb", header = F)
Data_post_20mb <- read.table("GATK_all_SNP_filtered_pass_cand_plusminus20mb_coverage_10kb", header = F)
Data_pre_20mb <- subset(Data_pre_20mb, V1 == "CM014754.1")
Data_post_20mb <- subset(Data_post_20mb, V1 == "CM014754.1")
Data_pre_20mb <- subset(Data_pre_20mb, V2 > 21226340 & V2 < 21841008)
Data_post_20mb <- subset(Data_post_20mb, V2 > 21226340 & V2 < 21841008)

Data_20mb <- rbind(Data_pre_20mb, Data_post_20mb)
Data_20mb$set <- c(rep("pre",62),rep("post",62))
Data_20mb$set <- factor(Data_20mb$set, levels = c("pre", "post"))
p4 <- ggplot(Data_20mb, aes(x=V2, y=V4, fill=set)) +  geom_bar(stat="identity", position=position_dodge(), show.legend = FALSE) +
  scale_fill_manual(values=c("lightgrey","darkblue")) + theme(aspect.ratio = 0.1) +
  theme_bw()

pdf("Coverage_stats.pdf", paper="a4", useDingbats=FALSE)
grid.arrange(p0,p1,p2,p3,p4,
             ncol=1)
dev.off()
