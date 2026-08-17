## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, March 2026.
## The goals of this R script are:
## To extract genomic outliers based on Tajima's D.

## 1. Setting up ----

## Packages:
library(tidyverse)

## Enforce no use of scientific number notation given large numbers of SNPs or chromosome windows:
options(scipen = 999)

## Working directory:
setwd("~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/TajimasD/0326/")

## Create directories to save outputs:
dir.create("data_outliers")

## 2. Function: Prepare scan data ----

## Define function:
get_TD <- function(taxon, pd) {

  ## Testing:
  #taxon <- "Pfi"
  
  ## Status:
  print(paste("Now processing", taxon, pd))
  
  ## Define window size:
  w <- 5000
  
  ## Tajima's D estimates:
  Db <- read.table(file = paste0("data_TajimasD/", taxon, "_B_w", w, ".Tajima.D"), header = TRUE)
  Dg <- read.table(file = paste0("data_TajimasD/", taxon, "_G_w", w, ".Tajima.D"), header = TRUE)
  
  ## In windows without SNPs to estimate Tajima's D, change not a number (NaN) to FALSE:
  Db <- Db[!is.nan(Db$TajimaD), ]
  Dg <- Dg[!is.nan(Dg$TajimaD), ]

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
  scan <- merge(Db, Dg, by = c("CHROM", "BIN_START", "BIN_END"), all = TRUE)
  
  ## Windows as numeric:
  scan$BIN_START <- as.numeric(scan$BIN_START)
  scan$BIN_END <- as.numeric(scan$BIN_END)
  
  ## Arrange:
  scan <- arrange(scan, CHROM, BIN_END)
  
  ## Get rid of windows with way too few SNPs, which can lead to unstable Tajima's D estimates:
  scan <- subset(scan, nSNPs_TDb >= 5 | nSNPs_TDg >= 5)
  
  ## Subtract Tajima's D green from brown:
  scan$deltaD <- scan$TDb - scan$TDg
  
  ## Absolute TD difference:
  scan$abs_deltaD <- abs(scan$deltaD)

  ## Define cutoffs and threshold (absolute delta TD):
	TD_topcut <- quantile(scan$abs_deltaD, pd, na.rm = TRUE)  ## Top cutoff.
  scan$outlier_TajimasD <- scan$abs_deltaD > TD_topcut ## ## Apply outlier cutoffs, top (absolute delta TD).
  
  ## Fix decimals to four digits:
  pd <- formatC(pd, format = "f", digits = 4)
  
  ## Outlier windows:
  outs <- scan[scan$outlier_TajimasD, ] ; dim(outs)
  
  ## BED format:
  bed <- outs[c("CHROM", "BIN_START", "BIN_END")]
  
  ## Convert window start to BED format (0-based start, 1-based end):
  bed$BIN_START <- bed$BIN_START - 1
  
  ## Remove NAs:
  bed <- bed[!is.na(bed$CHROM), ]
  
  ## Chromosome names:
  bed$CHROM <- factor(x = bed$CHROM, 
                      levels = sort(unique(bed$CHROM)),
                      labels = c(1:18, "Z"))
  
  ## Write BED file:
  write.table(bed,
              file = paste0("data_outliers/", taxon, "_outlier_windows_TajimasD_w", w, "_", pd, ".bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

} ## End of function.

## Outlier cutoff:
pd <- 0.9925

## Run:
for (t in c("Pcr", "Pfi", "Pga", "Pli", "Ppi", "Pva")) {
  get_TD(t, pd)
}

## End of script.
