#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### Lund University, Sweden.
### July 2024.
### The goals of this script are:
### To create and submit multiple SNP calling analyses, separately by chromosome and sample combination.

## Focal taxon:
#sp=Pcr
#sp=Pga
#sp=Pfi
#sp=Ppi
#sp=Pli
#sp=Pva

## Folder with raw sequence data:
## We'll use this to extract a list of specimens from the focal taxon.
#rawdir=/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/00_rawdata

## List samples of focal taxon in that folder:
#sample_list=$(ls ${rawdir} | grep ${sp})

## Loop over selected chromosomes:

## Allocate 4 hours:
#time=04:00:00
#for chrom in NC_011607.1 NC_041330.1 NC_041329.1 NC_041328.1 NC_041327.1

## Allocate 6 hours:
#time=06:00:00
#for chrom in NC_041326.1 NC_041325.1 NC_041324.1 NC_041323.1 NC_041322.1

## Allocate 7 hours:
#time=07:00:00
#for chrom in NC_041321.1 NC_041320.1 NC_041319.1 NC_041318.1 NC_041317.1

## Allocate 9 hours:
#time=09:00:00
#for chrom in NC_041316.1 NC_041315.1 NC_041314.1 NC_041313.1 NC_041312.1

## Redos:
time=04:00:00
for chrom in NC_011607.1

## Do:
do

## Status:
echo "Now doing ${chrom}!"

## Loop over listed samples to create and submit jobs:
#for sample in ${sample_list}

## Or, loop over selected samples:
for sample in \
Pli_GOD_H05

## Do:
do

## Status:
echo "Now doing ${sample}!"

## Create temporary job file from template:
cp \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling.job \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Replace target sample in temporary job file:
sed -i s/sample_field/${sample}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Replace target chromosome in temporary job file:
sed -i s/chrom_field/${chrom}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Replace target time in temporary job file:
sed -i s/time_field/${time}/ \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Submit job:
sbatch \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Remove temporary job file:
rm \
/cfs/klemming/projects/snic/naiss2024-23-107/Podarcis/scripts/03_SNP_calling_${sample}_${chrom}.tjob

## Close loops:
done
done

## Check:
squeue -u ivanpra | nl

## End of script.

