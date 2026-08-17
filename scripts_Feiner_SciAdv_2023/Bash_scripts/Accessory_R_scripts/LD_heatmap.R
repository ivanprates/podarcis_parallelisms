### This code creates fig. S6b and c

#dev.off()
rm(list=ls())
library(ggplot2)
library(tidyr)
library(gridExtra)
library(plotLD)
library(data.table)
library(LDheatmap)
library(tibble)
library(hrbrthemes)
library(dplyr)

Data <- read.table("1mb_pairwise_ld.ld", header = T)
Data <- Data[,-c(1,3,4,6)]

s <- min(Data[,1])
x <- 50000 #window size, can be changed (either 10000 or 50000)
k=1
Data_sum <- data.frame(matrix(ncol=4,nrow=10001, dimnames=list(NULL, c("midpoint_A", "midpoint_B", "R2_mean", "R2_95CI"))))

for(i in 1:50){ #number of blocks, can be changed (either 100 or 50)
  for(j in 1:50){ #number of blocks, can be changed (either 100 or 50)
    if(j>=i){
block <- Data[Data$BP_A>=s+(i-1)*x & Data$BP_A<(s+i*x) & Data$BP_B>=(s+j*x) & Data$BP_B<(s+(j+1)*x), ]
Data_sum$midpoint_A[k] <- (s + (i-1)*x + s + i*x)/2
Data_sum$midpoint_B[k] <- (s + j*x + s + (j+1)*x)/2
Data_sum$R2_mean[k] <- mean(block$R2)
Data_sum$R2_95CI[k] <- as.numeric(quantile(block$R2, 0.95))
j <- j+1
k <- k+1
  }}
i <- i+1  
}

Data_sum_noNA <- Data_sum[complete.cases(Data_sum), ]

pdf("1mb_pairwise_ld_50kb_blocks.pdf", width=5, height=5, useDingbats=FALSE) #change name!
ggplot(Data_sum_noNA,aes(x=midpoint_A/1000000,y=midpoint_B/1000000)) + theme_classic() +
  geom_tile(aes(fill=R2_95CI))+
  scale_fill_gradientn(colours=c("white","lightgrey","grey80","turquoise","deepskyblue3","blue3","navyblue","black"), limits=c(0,1), values=c(0,0.125,0.25,0.375,0.5,0.675,0.75,1), name="R2")+ xlab("position (Mb)") + ylab("position (Mb)") +
  scale_x_continuous(expand=c(0.02,0))+
  scale_y_continuous(expand=c(0.02,0)) +
  coord_fixed(ratio = 1)
dev.off()


###Plot the same, but with only the 60 individuals from IT
rm(list=ls())
Data <- read.table("1mb_pairwise_ld_60ind.ld", header = T)
Data <- Data[,-c(1,3,4,6)]

### plot windows
s <- min(Data[,1])
x <- 10000 #window size, can be changed (either 10000 or 50000)
k=1
Data_sum <- data.frame(matrix(ncol=4,nrow=10001, dimnames=list(NULL, c("midpoint_A", "midpoint_B", "R2_mean", "R2_95CI"))))

for(i in 1:100){ #number of blocks, can be changed (either 100 or 50)
  for(j in 1:100){ #number of blocks, can be changed (either 100 or 50)
    if(j>=i){
      block <- Data[Data$BP_A>=s+(i-1)*x & Data$BP_A<(s+i*x) & Data$BP_B>=(s+j*x) & Data$BP_B<(s+(j+1)*x), ]
      Data_sum$midpoint_A[k] <- (s + (i-1)*x + s + i*x)/2
      Data_sum$midpoint_B[k] <- (s + j*x + s + (j+1)*x)/2
      Data_sum$R2_mean[k] <- mean(block$R2)
      Data_sum$R2_95CI[k] <- as.numeric(quantile(block$R2, 0.95))
      j <- j+1
      k <- k+1
    }}
  i <- i+1  
}

Data_sum_noNA <- Data_sum[complete.cases(Data_sum), ]

pdf("1mb_pairwise_ld_60ind_10kb_blocks.pdf", width=5, height=5, useDingbats=FALSE) #change name!
ggplot(Data_sum_noNA,aes(x=midpoint_A/1000000,y=midpoint_B/1000000)) + theme_classic() +
  geom_tile(aes(fill=R2_95CI))+
  scale_fill_gradientn(colours=c("white","lightgrey","grey80","turquoise","deepskyblue3","blue3","navyblue","black"), limits=c(0,1), values=c(0,0.125,0.25,0.375,0.5,0.675,0.75,1), name="R2")+ xlab("position (Mb)") + ylab("position (Mb)") +
  scale_x_continuous(expand=c(0.02,0))+
  scale_y_continuous(expand=c(0.02,0)) +
  coord_fixed(ratio = 1)
dev.off()
