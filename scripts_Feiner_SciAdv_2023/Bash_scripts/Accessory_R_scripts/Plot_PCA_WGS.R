### This code creates fig. S1e-g

dev.off()
rm(list=ls())

library(ggplot2)
library(gridExtra)

#load all necessary input data files
Data_149ind <- read.csv("WGS_149ind_assignment.csv", sep=",", header=T)
evec <- read.table("populations_PCA.eigenvec") 
eval <- read.table("populations_PCA.eigenval")
evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
pve <- data.frame(PC = 1:20, pve = eval/sum(eval)*100)
#plot(pve$PC, pve$V1) #plot scree plot

#merge datasets
Data_149ind_PCs <- merge(Data_149ind, evec, by.x = "ID", by.y = "V1", all.x=T)

#Plot PCAs:

p1 <- ggplot(Data_149ind_PCs, aes(x=V3, y=V4, color=Lineage)) +
  geom_point(size=2) + scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

p2 <- ggplot(Data_149ind_PCs, aes(x=V3, y=V4, color=Greenness)) +
  geom_point(size=2) + scale_color_gradientn(colors = c("#816830","#675326","#4d3e1d","#336600","#339900","#66CC33")) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

p3 <- ggplot(Data_149ind_PCs, aes(x=V3, y=V4, color=Phenotype_group)) +
  geom_point(size=2) + scale_color_manual(values=c("#816830", "grey", "#66CC33")) +
  xlab(paste("PC1\n",evec1.pc, "% of observed genetic variation", sep="")) +
  ylab(paste("PC2\n",evec2.pc, "% of observed genetic variation", sep="")) +
  theme_bw() + theme(legend.position="none")

pdf("PCAs_149ind.pdf", width = 15  , height = 4.5, useDingbats=FALSE)
grid.arrange(p1,p2,p3, ncol=3)
dev.off()
