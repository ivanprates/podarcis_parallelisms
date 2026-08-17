#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### Lund University, Sweden.
### July 2024.
### The goals of this script are:
### To create and submit job files, each merging fastq files for those samples with more than two sequenced lanes.

## Looping over samples with more than one sequenced lane:
for sample in \
Pfi_CC_G105 \
Pfi_CC_G27 \
Pfi_CC_H46 \
Pfi_CI_A29 \
Pfi_CI_A82 \
Pfi_CO_A15 \
Pfi_CO_E02 \
Pfi_CO_G25

## Do:
do

## Status:
echo "Now submitting ${sample}!"

## Create temporary job file:
cp \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/00_merge_lanes.job \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/00_merge_lanes_${sample}.tjob

## Replace target control file in temporary job file:
sed -i s/sample_field/${sample}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/00_merge_lanes_${sample}.tjob

## Submit job:
sbatch \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/00_merge_lanes_${sample}.tjob

## Remove temporary job file:
rm \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/00_merge_lanes_${sample}.tjob

## Close loop:
done

## Check:
squeue -u ivanpra | nl

## End of script.

