#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### Lund University, Sweden.
### July 2024.
### The goals of this script are:
### To create and submit multiple analyses to process GVCF files, separately by taxon and chromosome.

## Focal taxon:
#sp=Pcr
sp=Pfi
#sp=Pga
#sp=Pli
#sp=Ppi
#sp=Pva

## Status:
echo "Now doing ${sp}!"

## Loop over selected chromosomes:
for chrom in \
NC_041326.1 \
NC_041325.1 \
NC_041324.1 \
NC_041323.1 \
NC_041322.1 \
NC_041321.1 \
NC_041320.1 \
NC_041319.1 \
NC_041318.1 \
NC_041317.1 \
NC_041316.1 \
NC_041315.1 \
NC_041314.1 \
NC_041313.1 \
NC_041330.1 \
NC_011607.1 \
NC_041327.1 \
NC_041328.1 \
NC_041329.1 \
NC_041312.1

## Do:
do

## Status:
echo "Now doing ${chrom}!"

## Create temporary job file from template:
cp \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs.job \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs_${sp}_${chrom}.tjob

## Replace target species in temporary job file:
sed -i s/sp_field/${sp}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs_${sp}_${chrom}.tjob

## Replace target chromosome in temporary job file:
sed -i s/chrom_field/${chrom}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs_${sp}_${chrom}.tjob

## Submit job:
sbatch \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs_${sp}_${chrom}.tjob

## Remove temporary job file:
rm \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/04_combine_GVCFs_${sp}_${chrom}.tjob

## Close loops:
done

## Check:
squeue -u ivanpra | nl

## End of script.

