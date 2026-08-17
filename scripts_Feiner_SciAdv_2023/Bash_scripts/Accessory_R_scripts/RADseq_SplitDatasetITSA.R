### This code also generates fig. S1 and S2

dev.off()
rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(openxlsx)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
options(scipen=999)
library(gridExtra)
library(conStruct)

#load all necessary input data files
Plink_fam <- read.csv("populations.plink_pruned.fam", sep=" ", header=F)
Phenotype_data <- read.csv("Phenotype_data_final_June2023.csv")
Pop_locations <- read.xlsx("PopsLocs.csv")
meanQ <- read.table("structure_result.2.meanQ", header=F, col.names = c("Q1","Q2"))
evec <- read.table("populations_PCA.eigenvec")
eval <- read.table("populations_PCA.eigenval")
evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
pve <- data.frame(PC = 1:20, pve = eval/sum(eval)*100)

#add phenotype data to .fam file and export (overwrites old .fam file)
Plink_fam_new <- merge(Plink_fam, Phenotype_data, by.x = "V2", by.y = "label", all.x=T, sort=F)
Plink_fam_new <- Plink_fam_new[order(Plink_fam_new$V2), ]
#add Q-values to phenptype data
Data <- cbind(Plink_fam_new,meanQ)
#add PC scores
Data <- merge(Data, evec, by.x = "V2", by.y = "V2", all.x=T, sort=F)
write.table(Data[,c(2,1,3:5,16)], file = "populations.plink_pruned_new.fam", sep = " ", col.names=F, quote=F, row.names=F)

#structure plot used in Extended Data Fig. 1a
Data_Q <- Data[,c(1,2,22,23)]
Q_sum <- as.data.frame(Data_Q %>% group_by(V1.x)  %>% dplyr::summarize(mean(Q2)))
colnames(Q_sum)[2] <- "Q2_mean"
Q_sum$Q1_mean <- 1-Q_sum$Q2_mean
Q_sum <- Q_sum[order(-Q_sum$Q2_mean), ]
Q_sum_long <- gather(Q_sum, region, Q_val, Q1_mean:Q2_mean, factor_key=TRUE)
pdf("Structure_plot.pdf", width=30, height=3, useDingbats = F)
ggplot(Q_sum_long, aes(x=V1.x, y=Q_val, fill=region)) +
  geom_bar(position="fill", stat="identity") +  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  scale_x_discrete(limits=Q_sum_long$V1.x[1:71]) + theme_bw() + theme(legend.position="none",
  panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#assign pops to lineage (IT and SA or mixed)

Data <- Data[,-c(7:9,11,12,14,15,17,18)]
#collapse to population-level data
data_summary <- as.data.frame(Data %>% group_by(V1.x)  %>% dplyr::summarize(mean(Q1)))

data_summary$V1.x[which(!data_summary$V1.x %in% Pop_locations$Abbreviation)] #None is missing
Data_sum <- merge(data_summary, Pop_locations, by.x="V1.x", by.y="Abbreviation", all=F, all.x=F, all.y=F)
colnames(Data_sum)[2] <- "meanQ"
Data_sum$Latitude <- as.numeric(Data_sum$Latitude)
Data_sum$Longitude <- as.numeric(Data_sum$Longitude)
Data_sum$region <- ifelse(Data_sum$meanQ >0.95, "IT",
                    ifelse(Data_sum$meanQ <0.05, "SA", "mixed"))
Data_sum %>% group_by(region) %>% tally()
write.table(Data_sum[,c(1:3,8:11,13)], file = "Pops.LineagesByQ.csv", sep = ",", col.names=T, quote=F, row.names=F)

Data_all <- merge(Data, Data_sum, by.x="V1.x", by.y="V1.x")

#this code generates Extended Data Fig. 1b and 1c
p1 <- ggplot(Data_all, aes(x=V3.y, y=V4.y, color=region)) +
  geom_point() + scale_color_manual(values=c("#E69F00","#D55E00", "#56B4E9")) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

p2 <- ggplot(Data_all, aes(x=V3.y, y=V4.y, color=V1.x)) +
  geom_point(alpha=0.2) + stat_ellipse() + scale_color_manual(values=rep(c("black"),71)) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

p3 <- ggplot(Data_all, aes(x=V3.y, y=V4.y, color=GreenResolved)) + 
  geom_point(alpha=0.8) + scale_color_gradientn(colors = c("#816830","#675326","#4d3e1d","#336600","#339900","#66CC33")) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

pdf("PCAs_890ind.pdf", width = 15  , height = 4.5, useDingbats=FALSE)
grid.arrange(p1,p2,p3, ncol=3)
dev.off()

#extract pure and mixed populations:
#SA <- All_samples_summ_with_location$V1[which(All_samples_summ_with_location$meanQ>0.9)] #SA pops
#IT <- All_samples_summ_with_location$V1[which(All_samples_summ_with_location$meanQ<0.1)] #IT pops
#All_samples_summ_with_location$V2[which((All_samples_summ_with_location$meanQ>0.1)&(All_samples_summ_with_location$meanQ<0.9))]

Data_all$sex <- gsub("F", "0", Data_all$sex)
Data_all$sex <- gsub("M", "1", Data_all$sex)

Data_all <- Data_all[complete.cases(Data_all[9]), ]

#export MyCov for common GWAS
#write.table(Data_all[,c(1,2,7,14:20)], file = "MyCov", sep = " ", col.names=F, quote=F, row.names=F)

Data_IT <- subset(Data_all, region =="IT")
Data_SA <- subset(Data_all, region =="SA")
Data_mixed <- subset(Data_all, region =="mixed")

#write.table(Data_all[,c(1,2)], file = "866ind", sep = " ", col.names=F, quote=F, row.names=F)
write.table(Data_IT[,c(1,2)], file = "ITind", sep = " ", col.names=F, quote=F, row.names=F)
write.table(Data_SA[,c(1,2)], file = "SAind", sep = " ", col.names=F, quote=F, row.names=F)
write.table(Data_mixed[,c(1,2)], file = "mixedind", sep = " ", col.names=F, quote=F, row.names=F)

### then run plink PCA on this subset of individuals
### see script_regionPC
###(plink -bfile populations.plink_pruned --keep ITind --pca --out IT/out --aec --make-bed)

install.packages("ggbiplot")
library(factoextra)

#for PC of 5 traits
Data_IT_sub <- Data_IT[complete.cases(Data_IT[,c(8:12)]), ]
IT_PC <- prcomp(Data_IT_sub[,c(8:12)], center = TRUE, scale. = TRUE)
print(IT_PC)
summary(IT_PC)
p3 <- fviz_pca_biplot(IT_PC, label = "var",col.var = "black", col.ind = "gray60", repel = T, title = "IT lineage")
Data_IT_sub <- cbind(Data_IT_sub,IT_PC$x[,1])
Data_IT_sub <- Data_IT_sub[,-c(1,3:47)]
Data_IT_PC <- merge(Data_IT, Data_IT_sub, by="V2",all.x=T)
write.table(Data_IT_PC[,c(2,1,3:5,48)], file = "IT_PC/out.fam", sep = " ", col.names=F, quote=F, row.names=F)

Data_SA_sub <- Data_SA[complete.cases(Data_SA[,c(8:12)]), ]
SA_PC <- prcomp(Data_SA_sub[,c(8:12)], center = TRUE, scale. = TRUE)
print(SA_PC)
summary(SA_PC)
p4 <- fviz_pca_biplot(SA_PC, label = "var",col.var = "black", col.ind = "gray60", repel = T, title = "SA lineage")
Data_SA_sub <- cbind(Data_SA_sub,SA_PC$x[,1])
Data_SA_sub <- Data_SA_sub[,-c(1,3:47)]
Data_SA_PC <- merge(Data_SA, Data_SA_sub, by="V2",all.x=T)
write.table(Data_SA_PC[,c(2,1,3:5,48)], file = "SA_PC/out.fam", sep = " ", col.names=F, quote=F, row.names=F)

Data_mixed_sub <- Data_mixed[complete.cases(Data_mixed[,c(8:12)]), ]
mixed_PC <- prcomp(Data_mixed_sub[,c(8:12)], center = TRUE, scale. = TRUE)
print(mixed_PC)
summary(mixed_PC)
Data_mixed_sub <- cbind(Data_mixed_sub,mixed_PC$x[,1])
Data_mixed_sub <- Data_mixed_sub[,-c(1,3:47)]
Data_mixed_PC <- merge(Data_mixed, Data_mixed_sub, by="V2",all.x=T)
write.table(Data_mixed_PC[,c(2,1,3:5,48)], file = "mixed_PC/out.fam", sep = " ", col.names=F, quote=F, row.names=F)

#this code creates Extended Data Fig. 2 
pdf("PCAs_890ind_biplot.pdf", width = 10, height = 8, useDingbats=FALSE)
grid.arrange(p3,p4, ncol=2)
dev.off()

#then read in PCA scores and make new Mycov file
evec_IT <- read.table("IT/out.eigenvec")
Data_IT <- Data_IT[,-c(13:44)]
Data_IT <- merge(Data_IT, evec_IT, by.x = "V2", by.y = "V2", sort=F)
write.table(Data_IT[,c(2,1,7,14:33)], file = "IT/MyCov", sep = " ", col.names=F, quote=F, row.names=F) #alt 2,1,7,14:33

evec_SA <- read.table("SA/out.eigenvec")
Data_SA <- Data_SA[,-c(13:44)]
Data_SA <- merge(Data_SA, evec_SA, by.x = "V2", by.y = "V2", sort=F)
write.table(Data_SA[,c(2,1,7,14:33)], file = "SA/MyCov", sep = " ", col.names=F, quote=F, row.names=F)

evec_mixed <- read.table("mixed/out.eigenvec")
Data_mixed <- Data_mixed[,-c(13:44)]
Data_mixed <- merge(Data_mixed, evec_mixed, by.x = "V2", by.y = "V2", sort=F)
write.table(Data_mixed[,c(2,1,17:33)], file = "mixed/MyCov", sep = " ", col.names=F, quote=F, row.names=F)
