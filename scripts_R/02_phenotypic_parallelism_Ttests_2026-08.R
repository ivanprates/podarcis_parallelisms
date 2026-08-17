#############################
## Script by Nathalie Feiner.

## Clean:
rm(list=ls())

## Packages:
library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(cowplot)
library(tidyverse)
library(grid)

## Directory:
setwd(dir = "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/")

## Trait data:
dp <- read.csv("morphological_analyses/data/Podarcis_greenbrowns_selected_traits_2026-03-15.csv")

## Traits to keep for downstream analyses:
selected <- c(
  "unique_identifier", "species", "group",
  "dorsum_total_color_h.theta",
  "dorsum_total_melanin_invlum",
  "right_OVS_UV_u",
  "right_OVS_UV_area_tot",
  "venter_total_melanin_area_rel",
  "svl",
  "headlength_by_SVL")

## Keep only these columns:
dp <- dp[, colnames(dp) %in% selected]
dp$svl <- as.numeric(dp$svl)

## Calculate angle (0 to 180):
angle <- function(a, b){
  return(acos(sum(a * b) / (sqrt(sum(a * a)) * sqrt(sum(b * b)) ) ) * (180/pi))
}

## Find all phenotype and species:
phenotype <- unique(dp$group)
species <- unique(dp$species)

## Total number of ecomorphs, species and combinations:
n.ecom <- length(phenotype)
n.isla <- length(species)

## Find all comparisons between ecomorphs:
phecomp <- t(combn(phenotype, 2))

## Total number of comparisons:
n.phecomp <- dim(phecomp)[1]

## Find all comparisions between species:
specomp <- t(combn(species, 2))

## total number of comparisons:
n.specomp <- dim(specomp)[1]

## built matrix for results:
theta <- deltaL <- array(NA, dim = c(n.phecomp,n.specomp))
colnames(theta) <- colnames(deltaL) <- apply(specomp, 1, paste, collapse = "-")
rownames(theta) <- rownames(deltaL) <- apply(phecomp, 1, paste, collapse = "-")

## Populate matrix:
for(i in 1:n.specomp){
  for(j in 1:n.phecomp){
    set1 <- subset(dp, group == phecomp[j, 1] & species %in% specomp[i, 1])
    set2 <- subset(dp, group == phecomp[j, 2] & species %in% specomp[i, 1])
    set3 <- subset(dp, group == phecomp[j, 1] & species %in% specomp[i, 2])
    set4 <- subset(dp, group == phecomp[j, 2] & species %in% specomp[i, 2])
    Spe1 <- Spe2 <- array(NA, dim = 7) ## Number of traits.
    for (k in 4:10) {
      if(length(set1[, k]) > 0 & length(set2[, k]) > 0){
        Spe1[k-3] <- as.numeric(t.test(set1[, k], set2[, k])$statistic)
      }else{
        Spe1[k-3] <- NA
      }
      if(length(set3[, k]) > 0 & length(set4[, k]) > 0){
        Spe2[k-3] <- as.numeric(t.test(set3[, k], set4[, k])$statistic)
      }else{
        Spe2[k-3] <- NA
      }
    }
    theta[j, i] <- angle(Spe1, Spe2) 
    deltaL[j, i] <- sqrt(sum(Spe1^2)) - sqrt(sum(Spe2^2)) ## Substract vector lengths in pairs.
  }
}
print("observed angles")
theta
print("observed deltaL")
deltaL

#### SIMULATIONS to test for parallel evolution:
n.sim <- 1000
perm_theta <- perm_deltaL <- array(NA, dim = c(n.phecomp, n.specomp, n.sim))
colnames(perm_theta) <- colnames(perm_deltaL) <- apply(specomp, 1, paste, collapse = "-")
rownames(perm_theta) <- rownames(perm_deltaL) <- apply(phecomp, 1, paste, collapse = "-")
perm_theta[, , 1] <- theta
perm_deltaL[, , 1] <- deltaL

for(i in 1:n.specomp){
  for(j in 1:n.phecomp){
    set1 <- subset(dp, group == phecomp[j, 1] & species %in% specomp[i, ])
    set2 <- subset(dp, group == phecomp[j, 2] & species %in% specomp[i, ])
    for(l in 2:n.sim){
      
      ## Permute the species:
      set1_random <- transform(set1, species = sample(species))
      set2_random <- transform(set2, species = sample(species))
      
      ## Split:
      set1Spe1 <- subset(set1_random, species %in% specomp[i, 1])
      set1Spe2 <- subset(set1_random, species %in% specomp[i, 2])
      set2Spe1 <- subset(set2_random, species %in% specomp[i, 1])
      set2Spe2 <- subset(set2_random, species %in% specomp[i, 2])
      Spe1 <- array(NA, dim = 7)
      Spe2 <- array(NA, dim = 7)
      for (k in 4:10) {
        if(length(set1Spe1[, k]) > 1 & length(set2Spe1[, k]) > 1){
          Spe1[k-3] <- as.numeric(t.test(set1Spe1[, k], set2Spe1[, k])$statistic)
        } else {
          Spe1[k-3] <- NA
        }
        if(length(set1Spe2[, k]) > 1 & length(set2Spe2[, k]) > 1){
          Spe2[k-3] <- as.numeric(t.test(set1Spe2[, k], set2Spe2[, k])$statistic)
        } else {
          Spe2[k-3] <- NA
        }
      }
      perm_theta[j, i, l] <- angle(Spe1, Spe2) 
      perm_deltaL[j, i, l] <- sqrt(sum(Spe1^2)) - sqrt(sum(Spe2^2)) ## Subtract vector lengths in pairs.
    }
  }
}

## Upper boundary:
perm_theta_975CI <- apply(perm_theta, c(1, 2), quantile, na.rm = TRUE, probs = 0.975)
perm_deltaL_975CI <- apply(perm_deltaL, c(1, 2), quantile, na.rm = TRUE, probs = 0.975)
perm_deltaL_025CI <- apply(perm_deltaL, c(1, 2), quantile, na.rm = TRUE, probs = 0.025)
perm_deltaL_mean <- apply(perm_deltaL, c(1, 2), mean, na.rm = TRUE)

## Test if observed data is consistent with parallelism:
parallel_theta <- perm_theta_975CI - theta
parallel_deltaL <- perm_deltaL_975CI - deltaL
print("parallel_theta")
parallel_theta
print("sig parallel theta?")
theta <= perm_theta_975CI
print("parallel_deltaL")
parallel_deltaL
print("sig parallel deltaL?")
deltaL >= perm_deltaL_025CI & deltaL <= perm_deltaL_975CI

## Get P-values:
P_theta <- P_deltaL <- array(NA, dim = c(n.phecomp, n.specomp))
colnames(P_theta) <- colnames(P_deltaL) <- apply(specomp, 1, paste, collapse = "-")
rownames(P_theta) <- rownames(P_deltaL) <- apply(phecomp, 1, paste, collapse = "-")

for(i in 1:n.specomp){
  for(j in 1:n.phecomp){
    P_theta[j,i] <- 1-sum(as.vector(na.omit(perm_theta[j,i,2:n.sim])) < theta[j,i])/length(as.vector(na.omit(perm_theta[j,i,2:n.sim]))) #is obs angle different from parallel?
    if(deltaL[j,i]>perm_deltaL_mean[j,i] & !is.na(deltaL[j,i])) 
      P_deltaL[j,i] <- 1-sum(as.vector(na.omit(perm_deltaL[j,i,2:n.sim])) < deltaL[j,i])/length(as.vector(na.omit(perm_deltaL[j,i,2:n.sim])))
    else
      P_deltaL[j,i] <- 1-sum(as.vector(na.omit(perm_deltaL[j,i,2:n.sim])) > deltaL[j,i])/length(as.vector(na.omit(perm_deltaL[j,i,2:n.sim])))
  }
}

print("P-value theta")
P_theta
print("P-value deltaL")
P_deltaL

#### Test for orthogonal evolution:
n.sim <- 1000 
boot_theta <- array(NA, dim = c(n.phecomp, n.specomp, n.sim))
colnames(boot_theta) <- apply(specomp, 1, paste, collapse = "-")
rownames(boot_theta) <- apply(phecomp, 1, paste, collapse = "-")
boot_theta[,,1] <- theta

for(i in 1:n.specomp){
  for(j in 1:n.phecomp){
    phe1Spe1 <- subset(dp, group == phecomp[j,1] & species %in% specomp[i,1])
    phe2Spe1 <- subset(dp, group == phecomp[j,2] & species %in% specomp[i,1])
    phe1Spe2 <- subset(dp, group == phecomp[j,1] & species %in% specomp[i,2])
    phe2Spe2 <- subset(dp, group == phecomp[j,2] & species %in% specomp[i,2])
    for(l in 2:n.sim){
      #bootstrap data
      phe1Spe1_boot <- sample_frac(phe1Spe1, 1, replace=T)
      phe2Spe1_boot <- sample_frac(phe2Spe1, 1, replace=T)
      phe1Spe2_boot <- sample_frac(phe1Spe2, 1, replace=T)
      phe2Spe2_boot <- sample_frac(phe2Spe2, 1, replace=T)
      Spe1 <- array(NA, dim = 7)
      Spe2 <- array(NA, dim = 7)
      for (k in 4:10){
        if(length(phe1Spe1_boot[,k])>1 & length(phe2Spe1_boot[,k])>1 & sd(phe1Spe1_boot[,k])>0 & sd(phe2Spe1_boot[,k])>0) {
          Spe1[k-3] <- as.numeric(t.test(phe1Spe1_boot[,k],phe2Spe1_boot[,k])$statistic)
        }else{
          Spe1[k-3] <- NA
        }
        if(length(phe1Spe2_boot[,k])>1 & length(phe2Spe2_boot[,k])>1 & sd(phe1Spe2_boot[,k])>0 & sd(phe2Spe2_boot[,k])>0) {
          Spe2[k-3] <- as.numeric(t.test(phe1Spe2_boot[,k],phe2Spe2_boot[,k])$statistic)
        }else{
          Spe2[k-3] <- NA
        }
      }
      boot_theta[j,i,l] <- angle(Spe1,Spe2) 
    }
  }
}

## Upper boundary:
boot_theta_975CI <- apply(boot_theta, c(1,2), quantile, na.rm=TRUE, probs=0.975)

## Lower boundary:
boot_theta_025CI <- apply(boot_theta, c(1,2), quantile, na.rm=TRUE, probs=0.025)

## Test if bootstrapped data includes 90 degrees:
print("bootstrapping of theta - upper CI includes 90?")
boot_theta_975CI - 90 ## Negative values mean that it doesn't include 90 - reject orthogonal/random evolution.
print("Random angle?")
boot_theta_975CI >= 90

## Get P-values:
P_theta_boot <- array(NA, dim = c(n.phecomp, n.specomp))
colnames(P_theta_boot) <- apply(specomp, 1, paste, collapse = "-")
rownames(P_theta_boot) <- apply(phecomp, 1, paste, collapse = "-")

for(i in 1:n.specomp){
  for(j in 1:n.phecomp){
    P_theta_boot[j, i] <- 1-sum(as.vector(na.omit(boot_theta[j, i, 2:n.sim])) < 90, na.rm = TRUE)/length(as.vector(na.omit(boot_theta[j,i,2:n.sim]))) #is obs angle different from parallel?
  }
}

print("P-value theta bootstrapped (incl 90?)")
P_theta_boot

###Plotting data to check if patterns seem to make sense:
#plot the distribution of permuted angles simulating parallel evolution (green) and bootstrapped observed angle (red)
plots <- list()
k <- 1
for(i in 1:21){ 
  for(j in 1){
    if(!is.na(theta[j, i])){
      parallel <- as.data.frame(perm_theta[j, i, 1:n.sim])
      bootstrap <- as.data.frame(boot_theta[j, i, 1:n.sim])
      colnames(parallel)[1] <- colnames(bootstrap)[1] <- "angle"
      angles <- rbind(parallel, bootstrap)
      angles[1:1000, 2] <- "parallel"
      angles[1001:2000, 2] <- "bootstrap"
      colnames(angles)[2] <- "Dataset"
      
      plots[[k]] <- ggplot(angles, aes(x = angle, fill = Dataset)) +
        geom_histogram(alpha = 0.7, position = "identity", aes(y = ..density..), color = "grey") +
        geom_vline(xintercept = quantile(parallel[,1], c(0.975), na.rm = T), size = .4, colour="blue")+ 
        geom_vline(xintercept = quantile(parallel[,1], c(.025), na.rm=T), size=.4, colour="blue")+ 
        geom_vline(xintercept = quantile(bootstrap[,1], c(.975), na.rm=T), size=.4, colour="red")+ 
        geom_vline(xintercept = quantile(bootstrap[,1], c(.025), na.rm=T), size=.4, colour="red")+ 
        geom_vline(xintercept = theta[j,i], size=1, colour="black") +
        geom_density(alpha=0.7, size=0.3) + xlim(0, 90) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),plot.title = element_text(size = 10)) +
        theme(plot.margin=unit(c(0,0,0.1,0),"cm")) +
        labs(title = (paste(colnames(perm_theta)[i],rownames(perm_theta)[j]))) 
      k <- k + 1
    }
  }
} 

## Store:
#pdf("FinalDataset_angles.pdf", width=10, height=14) #change name!
#do.call(plot_grid, c(plots, ncol = 8))

## Move to a new page:
grid.newpage()

## Create layout:
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 3)))

## A helper function to define a region on the layout:
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}

## Arrange the plots:
print(plots[[1]], vp = define_region(row = 1, col = 1)) ## Span over 6 columns.
print(plots[[2]], vp = define_region(row = 2, col = 1))
print(plots[[3]], vp = define_region(row = 3, col = 1))
print(plots[[4]], vp = define_region(row = 4, col = 1))
print(plots[[5]], vp = define_region(row = 5, col = 1))
print(plots[[6]], vp = define_region(row = 6, col = 1))
print(plots[[7]], vp = define_region(row = 7, col = 1))
print(plots[[8]], vp = define_region(row = 1, col = 2))
print(plots[[9]], vp = define_region(row = 2, col = 2))
print(plots[[10]], vp = define_region(row = 3, col = 2))
print(plots[[11]], vp = define_region(row = 4, col = 2))
print(plots[[12]], vp = define_region(row = 5, col = 2))
print(plots[[13]], vp = define_region(row = 6, col = 2))
print(plots[[14]], vp = define_region(row = 7, col = 2))
print(plots[[15]], vp = define_region(row = 1, col = 3))
print(plots[[16]], vp = define_region(row = 2, col = 3))
print(plots[[17]], vp = define_region(row = 3, col = 3))
print(plots[[18]], vp = define_region(row = 4, col = 3))
print(plots[[19]], vp = define_region(row = 5, col = 3))
print(plots[[20]], vp = define_region(row = 6, col = 3))
print(plots[[21]], vp = define_region(row = 7, col = 3))
#dev.off()

## Summarize the results and plot pairwise comparisons:
res_theta_parallel <- theta <= perm_theta_975CI  ## Which comparisons are indistinguishable from parallel.
res_theta_random <- boot_theta_975CI >= 90 ## Which comparisons are indistinguishable from 90 degrees (i.e., random).
res_deltaL <- deltaL >= perm_deltaL_025CI & deltaL <= perm_deltaL_975CI ## Which comparisons have indistinguishable length.

summary_theta <- data.frame(
  comparison   = colnames(res_theta_parallel),
  angle_deg    = round(as.numeric(theta[1, ]), digits=1),
  parallel     = as.logical(res_theta_parallel[1, ]),
  orthogonal   = as.logical(res_theta_random[1, ]),
  same_length  = as.logical(res_deltaL[1, ])
) %>%
  mutate(
    angle = case_when(
      parallel ~ "parallel",
      orthogonal ~ "not parallel",
      TRUE ~ "moderately parallel"
    ),
    length = case_when(
      same_length ~ "same length",
      TRUE ~ "different length"
    )
  ) %>%
  select(comparison, angle_deg, angle, length)

summary_theta

## Split comparison into two species columns:
summary_theta <- summary_theta %>%
  separate(comparison, into = c("sp1", "sp2"), sep = "-", remove = FALSE)

## Define species order:
all_species <- c("muralis", "filfolensis", "pityusensis", "vaucheri", "liolepis", "cretensis", "gaigeae")
species_labels <- c("P. muralis", "P. filfolensis", "P. pityusensis", "P. vaucheri", "P. liolepis", "P. cretensis", "P. gaigeae")
names(species_labels) <- all_species

## Factor with consistent levels:
summary_theta$sp1 <- factor(summary_theta$sp1, levels = all_species)
summary_theta$sp2 <- factor(summary_theta$sp2, levels = all_species)

summary_theta_full <- bind_rows(
  summary_theta,
  summary_theta %>% rename(sp1 = sp2, sp2 = sp1)
)

summary_theta_lower <- summary_theta_full %>%
  filter(as.integer(factor(sp1, levels = all_species)) > 
           as.integer(factor(sp2, levels = all_species)))

#pdf("./plots/Parallelism_pheno_angles_V1.pdf", width = 6, height = 6) ## Change name!
ggplot(summary_theta_lower, aes(x = sp2, y = sp1, fill = angle)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(angle_deg, 0), "\u00B0")), size = 3.5, color = "black", fontface = "bold") +
  scale_fill_manual(
    values = c("parallel"            = "#0AD8E1",
               "moderately parallel" = "#86A8E7",
               "not parallel"        = "#D16B74"),
    name = NULL
  ) +
  scale_x_discrete(limits = all_species, labels = species_labels) +
  scale_y_discrete(limits = rev(all_species), labels = species_labels) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(face = "italic"),
    panel.grid  = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.8),
  )
#dev.off()

mean(summary_theta_lower$angle_deg) ## Extract mean angle.
sd(summary_theta_lower$angle_deg)   ## And its std.

## Save results into supplementary table:
write.csv(summary_theta_full, file = "drafts/Table_S2_parallelism_analyses.csv", row.names = FALSE)

## End of script.
