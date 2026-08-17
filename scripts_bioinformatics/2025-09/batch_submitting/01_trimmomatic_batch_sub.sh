#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### Lund University, Sweden.
### July 2024.
### The goals of this script are:
### To create and submit multiple Trimmomatic analyses.

## Focal taxon:
#sp=Pcr
#sp=Pga
sp=Pfi
#sp=Ppi
#sp=Pli
#sp=Pva

## Folder with raw sequence data:
## We'll use this to extract a list of specimens from the focal taxon.
rawdir=/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/00_rawdata

## List samples of focal taxon in that folder:
sample_list=$(ls ${rawdir} | grep ${sp})

## Loop over these samples to create and submit jobs:
for sample in ${sample_list}

## Or, loop over selected samples:
#for sample in \
#Pva_TZ_K90

## Do:
do

## Status:
echo "Now submitting ${sample}!"

## Create temporary job file from template:
cp \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/01_trimmomatic.job \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/01_trimmomatic_${sample}.tjob

## Replace sample in temporary job file:
sed -i s/sample_field/${sample}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/01_trimmomatic_${sample}.tjob

## Submit job:
sbatch \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/01_trimmomatic_${sample}.tjob

## Remove temporary job file:
rm \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/01_trimmomatic_${sample}.tjob

## Close loop:
done

## Check:
squeue -u ivanpra

## End of script.

