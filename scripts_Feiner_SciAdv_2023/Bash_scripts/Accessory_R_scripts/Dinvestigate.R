### This code creates fig. S3

library(ggplot2)
library(tidyr)
library(gridExtra)

# read in the results from Dinvestigate
Data <- read.table("P1_P2_P3_localFstats__100_50.txt",as.is=T,header=T)

p1 <- ggplot(Data, aes(x=windowStart, y=D)) +
  geom_point() + xlim(c(21200000,21900000)) +
  xlab(paste("Scaffold 12")) +
  ylab(paste("f_dM")) + geom_hline(yintercept=quantile(Data$D, 0.9975), linetype="dashed", color = "red") +
  geom_vline(xintercept=21342501) + geom_vline(xintercept=21370000) +
  theme_bw() + ggtitle("D")

p2 <- ggplot(Data, aes(x=windowStart, y=f_dM)) +
  geom_point() + xlim(c(21200000,21900000)) +
  xlab(paste("Scaffold 12")) +
  ylab(paste("f_dM")) + geom_hline(yintercept=quantile(Data$f_dM, 0.9975), linetype="dashed", color = "red") +
  geom_vline(xintercept=21342501) + geom_vline(xintercept=21370000) +
  theme_bw() + ggtitle("f_dM")

pdf("Dsuite_output.pdf", width=12, height=5, useDingbats=FALSE)
grid.arrange(p1,p2, ncol=1)
dev.off()
