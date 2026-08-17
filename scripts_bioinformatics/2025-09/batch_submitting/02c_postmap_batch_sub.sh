#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### Lund University, Sweden.
### July 2024.
### The goals of this script are:
### To create and submit multiple genome mapping analyses.

## Focal taxon:
#sp=Pcr
#sp=Pga
#sp=Ppi
#sp=Pli
#sp=Pva

## Folder with raw sequence data:
## We'll use this to extract a list of specimens from the focal taxon.
#rawdir=/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/00_rawdata

## List samples of focal taxon in that folder:
#sample_list=$(ls ${rawdir} | grep ${sp})

## List samples in taxon-abbpop combination in folder:
#sample_list=$(ls ${rawdir} | grep Pfi_CI)

## Loop over these samples to create and submit jobs:
#for sample in ${sample_list}

## Alternatively, loop over selected samples:
for sample in \
Pfi_CO_A15 \
Pfi_CO_A78 \
Pfi_CO_E02 \
Pfi_CO_E22 \
Pfi_CO_F73 \
Pfi_CO_G101 \
Pfi_CO_G105 \
Pfi_CO_G25 \
Pfi_CO_G27 \
Pfi_CO_G58 \
Pfi_CO_H103 \
Pfi_CO_H16 \
Pfi_CO_H46 \
Pfi_CO_H48 \
Pfi_CO_H71

## Do:
do

## Status:
echo "Now submitting ${sample}!"

## Create temporary job file from template:
cp \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/02c_postmap.job \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/02c_postmap_${sample}.tjob

## Replace target sample in temporary job file:
sed -i s/sample_field/${sample}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/02c_postmap_${sample}.tjob

## Submit job:
sbatch \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/02c_postmap_${sample}.tjob

## Remove temporary job file:
rm \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/02c_postmap_${sample}.tjob

## Close loop:
done

## Check:
squeue -u ivanpra | nl

## End of script.
