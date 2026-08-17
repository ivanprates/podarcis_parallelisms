rm(list=ls())

library(openxlsx)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(pcaMethods)
library(RRPP)
library(dplyr)
library(scales)

options(scipen=999)

#Import the current new data set - this includes data up till 2018
data <- read.csv("Input_files/Phenotype_data_final_June2023.csv")

#change some colnames
colnames(data)[c(11,14,15,16)] <- c("green","black","blue","rel.headlength")

############### Fig. 1E,F - Phenotypic trajectories along transect (SA - hybrid zone - IT)

#cline along the coast
clinecoast <- subset(data, abbpop == "LO" | abbpop == "NL" | abbpop == "VA" | 
            abbpop == "ME" | abbpop == "GN" | abbpop == "RA" | abbpop == "SL" | abbpop == "LE" | 
            abbpop == "ST" | abbpop == "MG" | abbpop == "VI" | abbpop == "CA" | abbpop == "CR" |
            abbpop == "CN" | abbpop == "VE") 
cline_loc <- data.frame(abbpop = c("LO","NL","VA","ME","GN","RA","SL","LE","ST",
                                   "MG","VI","CA","CR","CN","VE"), 
                        pop_order = c("A_LO","B_NL","C_VA","D_ME","E_GN","F_RA","G_SL",
                                      "H_LE","I_ST","J_MG","K_VI","L_CA","M_CR","N_CN","O_VE"))

clinecoast <- merge(cline_loc,clinecoast, by ="abbpop") 
attach(clinecoast)
clinecoast <- clinecoast[order(pop_order),]
detach(clinecoast)

clinecoast <- clinecoast[complete.cases(clinecoast[,c(9,17,12,15,16)]), ]
traits <- as.matrix(scale(clinecoast[,c(9,17,12,15,16)]))
#traits_col <- as.matrix(scale(clinecoast[,c(12,15,16)]))

fit <- lm.rrpp(traits ~ pop_order*sex, data = clinecoast, iter = 999)
reveal.model.designs(fit)
TA <- trajectory.analysis(fit, groups = clinecoast$sex,
                          traj.pts = clinecoast$pop_order, print.progress = F)
summary(TA, attribute = "MD", show.trajectories=TRUE) # Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "TC", angle.type = "deg", show.trajectories=TRUE) # Correlations (angles) between trajectories
summary(TA, attribute = "SD", show.trajectories=TRUE) # No shape differences between vectors

# Retain results
TA.summary <- summary(TA, attribute = "SD")
TA.summary$summary.table
# Plot results
pdf("Output_plots/Trajectory.pdf", height = 5, useDingbats=FALSE)
TP <- plot(TA, pch = c(3,5)[as.numeric(as.factor(clinecoast$sex))], 
           col = c("darkgrey", "lightgrey")[as.numeric(as.factor(clinecoast$sex))],
           cex = 0.8)
add.trajectories(TP, traj.cex=2.2, traj.pch = 21, start.bg = "#5BB3E4", end.bg="#E59F24")
dev.off()
#legend("topright", paste(c(unique(as.factor(clinecoast$abbpop)))), pch = 21)

summary.trajectory.analysis(TA, attribute ="MD",show.trajectories = TRUE)

# Predictions (holding alternative effects constant)
shapeDF <- data.frame(sex = c(rep("M",17),rep("F",17)), 
                              pop_order = c(rep(c("A_LO","B_NL","C_SO","D_SM","E_VA","F_ME","G_GN","H_RA","I_SL",
                                                 "J_LE","K_ST","L_MG","M_VI","N_CA","O_CR","P_CN","Q_VE"),2))) # levels(colourApp$sex), alt = 
shapePreds <- predict(fit, shapeDF, confidence = 0.95)
summary(shapePreds)
summary(shapePreds, PC = TRUE)
shapePreds$pca
# Plot prediction
plot(shapePreds, PC = TRUE)
plot(shapePreds, PC = TRUE, ellipse = TRUE)

plot(shapePreds, PC = TRUE, ellipse = TRUE,
     pch = 19, col = 1:17)


###population means
mphen<-aggregate(cbind(msvl=clinecoast$svl, mmass=clinecoast$mass, mhl=clinecoast$rel.headlength, 
                       mgreen=clinecoast$green, mblack=clinecoast$black, mblue=clinecoast$blue),
                 by=list(abbpop=clinecoast$pop_order, sex=clinecoast$sex), mean,na.rm=TRUE)
mphen$pop_order <- as.integer(as.factor(mphen$abbpop))
mphen_M <- subset(mphen, sex=="M")
mphen_F <- subset(mphen, sex=="F")

mphen_M$mgreen_01 <- rescale(mphen_M$mgreen)

ggplot(mphen_M, aes(x=pop_order, y=mgreen_01)) + geom_point() + 
  geom_smooth(method="glm", method.args=list(family="quasibinomial"), se=FALSE)
