## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, August 2026.
## The goals of this R script are:
## To compose manuscript supplementary figure S1.

## 1. Getting ready ----

## Packages:
library(ggplot2)
library(patchwork)
library(tidyverse)

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Path to working directory: Ivan's Linux partition:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/"

## Set working directory:
setwd(path)

## Outlier threshold:
cutoff <- 0.9925

## Loop over species:
for (sp in c("Pcr", "Pga", "Pfi", "Pli", "Ppi", "Pva")) {

## Species:
#sp <- "Pcr"
#sp <- "Pga"
#sp <- "Pfi"
#sp <- "Pli"
#sp <- "Ppi"
#sp <- "Pva"

## Species label for plot:
if (sp == "Pcr") { sp_label <- "Podarcis cretensis" }
if (sp == "Pga") { sp_label <- "Podarcis gaigeae" }
if (sp == "Pfi") { sp_label <- "Podarcis filfolensis" }
if (sp == "Pli") { sp_label <- "Podarcis liolepis" }
if (sp == "Ppi") { sp_label <- "Podarcis pityusensis" }
if (sp == "Pva") { sp_label <- "Podarcis vaucheri" }

## Status:
print(paste0("Now processing ", sp_label))

## 2. Plot Tajima's D ----

## Define window size:
w <- 5000

## Tajima's D estimates:
Db <- read.table(file = paste0("TajimasD/0326/data_TajimasD/", sp, "_B_w", w, ".Tajima.D"), header = TRUE)
Dg <- read.table(file = paste0("TajimasD/0326/data_TajimasD/", sp, "_G_w", w, ".Tajima.D"), header = TRUE)
  
## In windows without SNPs to estimate Tajima's D, change not a number (NaN) to NA:
Db$TajimaD[is.nan(Db$TajimaD)] <- NA
Dg$TajimaD[is.nan(Dg$TajimaD)] <- NA
  
## Column names:
names(Db)[names(Db) == "TajimaD"] <- "TDb"
names(Dg)[names(Dg) == "TajimaD"] <- "TDg"
names(Db)[names(Db) == "N_SNPS"] <- "nSNPs_TDb"
names(Dg)[names(Dg) == "N_SNPS"] <- "nSNPs_TDg"
  
## Bins for Tajima's D data:
Db$BIN_END <- Db$BIN_START + w
Dg$BIN_END <- Dg$BIN_START + w
Db$BIN_START <- Db$BIN_START + 1
Dg$BIN_START <- Dg$BIN_START + 1
  
## Merge D data:
scanTD <- merge(Db, Dg, by = c("CHROM", "BIN_START", "BIN_END"), all = TRUE)
  
## Windows as numeric:
scanTD$BIN_START <- as.numeric(scanTD$BIN_START)
scanTD$BIN_END <- as.numeric(scanTD$BIN_END)

## Define factor for ordered chromosomes, chromosome labels in plots:
scanTD$CHROM <- factor(
                  x = scanTD$CHROM, 
                  levels = sort(unique(scanTD$CHROM)),
                  labels = c(1:18, "Z"))

## Arrange:
scanTD <- arrange(scanTD, CHROM, BIN_END)
  
## Subtract Tajima's D green from brown:
scanTD$deltaD <- scanTD$TDb - scanTD$TDg

## Absolute TD difference:
scanTD$abs_deltaD <- abs(scanTD$deltaD)

## Mid position of windows:
scanTD$win_mid <- (( scanTD$BIN_START+ scanTD$BIN_END ) /2 )

## chromosome lengths:
chr_TD <- scanTD %>% 
  group_by(CHROM) %>% 
  summarise(
    chr_max = max(BIN_END), 
    chr_center = (( min(BIN_START) + max(BIN_END) ) / 2 ),
    .groups = "drop")

## Cumulative chrom start position:
chr_TD$cum_start <- cumsum(lag(chr_TD$chr_max, default = 0))

## Cumulative chrom mid position:
chr_TD$cum_mid <- chr_TD$cum_start + chr_TD$chr_center

## Add chrom length info to scanTD data:
scanTD <- left_join(scanTD, chr_TD, by = "CHROM")

## Cumulative window middle:
scanTD$cum_mid <- scanTD$win_mid + scanTD$cum_start

## Cutoff:
cutoff_TD <- quantile(scanTD$abs_deltaD, as.numeric(cutoff), na.rm = TRUE)

## Plot scan of Tajima' D values:
t <- ggplot() +
  
  ## Points:
  geom_point(data = scanTD, size = 0.75, shape = 16, aes(x = cum_mid, y = abs_deltaD, color = CHROM)) +
  scale_color_manual(values = rep(c("grey70", "grey30"), length.out = length(unique(scanTD$CHROM)))) +
  
  ## Cutoff line:
  geom_hline(yintercept = cutoff_TD, color = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.35) +
    
  ## Scales, labels:
  labs(y = "|Δ Tajima's D|", x = "Chromosome") +
  scale_x_continuous(
    limits = range(scanTD$cum_mid),
    breaks = chr_TD$cum_mid,
    labels = chr_TD$CHROM,
    expand = expansion(mult = c(0.01, 0.01))) +
  
  ## Appearance:
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.line.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10, colour = "black"),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0))
  
## 3. XtX scans ----

## Read XtX estimates:
scanXtX <- read.table(file = paste0("baypass/0326/data_XtX/", sp, "_median_summary_pi_xtx.out"), header = TRUE)
names(scanXtX)[names(scanXtX) == "M_XtX"] <- "XtX"

## SNP information:
SNPinfo <- read.table(file = paste0("baypass/0326/data_XtX/", sp, "_snp_ids.txt"), header = FALSE)
names(SNPinfo) <- c("CHROM", "SNP", "MRK")

## Combine:
scanXtX_SNPs <- merge(scanXtX, SNPinfo)
names(scanXtX_SNPs)

## Window size for plotting:
## Using windows for plotting to become computationally tractable.
win_size <- 1000

## Bin into windows:
scan_windows <- scanXtX_SNPs %>%
  group_by(CHROM) %>%
   mutate(
    window = cut(SNP,
                 breaks = seq(floor(min(SNP)/win_size)*win_size,
                              ceiling(max(SNP)/win_size)*win_size,
                              by = win_size),
                 labels = FALSE,
                 include.lowest = TRUE),
    win_start = floor(min(SNP)/win_size)*win_size + (window - 1) * win_size + 1,
    win_end = floor(min(SNP)/win_size)*win_size + window * win_size)

## Estimate average XtX per window, annotate whether window includes outlier SNPs:
scan_windows <- scan_windows %>%
group_by(CHROM, window, win_start, win_end) %>%
  summarise(
    XtX_mean = mean(XtX, na.rm = TRUE),
    .groups = "drop")

## Mid position of windows:
scan_windows$win_mid <- (( scan_windows$win_start + scan_windows$win_end) / 2 )

## chromosome lengths:
chr_XtX <- scan_windows %>% 
  group_by(CHROM) %>% 
  summarise(
    chr_max = max(win_end), 
    chr_center = (( min(win_start) + max(win_end) ) / 2 ),
    .groups = "drop")

## Cumulative chrom start position:
chr_XtX <- chr_XtX %>% arrange(CHROM)
chr_XtX$cum_start <- cumsum(lag(chr_XtX$chr_max, default = 0))

## Add chrom length info to scan data:
scan_windows <- left_join(scan_windows, chr_XtX, by = "CHROM")

## Cumulative window middle:
scan_windows$cum_mid <- scan_windows$win_mid + scan_windows$cum_start

## Confirm we have 19 CHROM:
length(unique(scan_windows$CHROM))

## Define factor for ordered chromosomes, chromosome labels in plots:
scan_windows$CHROM <- factor(x = scan_windows$CHROM, 
                             levels = sort(unique(scan_windows$CHROM)),
                             labels = c(1:18, "Z"))

## Window-based cutoff:
cutoff_XtX <- quantile(scan_windows$XtX_mean, as.numeric(cutoff), na.rm = TRUE)

## Baypass cutoff:
BP_cutoff_XtX <- read.table(paste0("baypass/0326/data_XtX/", sp, "_mean_xtx_threshold_", cutoff, ".txt"), header = FALSE)

## Plot scan of XtX values:
x <- ggplot() +
  
  ## Plot windowed XtX:
  geom_point(data = scan_windows, 
    size = 0.5, shape = 16,
    mapping = aes(x = cum_mid, y = XtX_mean, color = CHROM)) +
  scale_color_manual(values = rep(c("grey70", "grey30"), length.out = length(unique(scan_windows$CHROM)))) +
  
  ## Cutoff line from window distribution:
  #geom_hline(yintercept = cutoff_XtX, color = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.35) +
  
  ## Cutoff from BayPass:
  geom_hline(yintercept = BP_cutoff_XtX$V1, color = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.35) +
  
  ## Axis:
  scale_x_continuous(limits = range(scan_windows$cum_mid), expand = expansion(mult = c(0.01, 0.01))) +
  
  ## Appearance:
  labs(y = "XtX statistics", title = sp_label) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.line.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0))

## 4. Compose figure (Fig. S1) ----

## Top part:
FigS1 <- ( x / t )

## For species:
ggsave(plot = FigS1, width = 40, height = 10, units = "cm", dpi = 300,
       filename = paste0(path, "figures/FigS1_", sp, "_w", win_size, ".jpg"))

} ## Close loop (species).

## End of script.
