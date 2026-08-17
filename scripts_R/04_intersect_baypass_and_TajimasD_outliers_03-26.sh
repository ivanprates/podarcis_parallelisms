## R script by Ivan Prates (ivanprates.org) ----
## Lund University, Sweden, March 2026.
## The goals of this shell script are:
## To overlay BayPass SNPs and Tajima's D windows:
## To extract genes sitting on overlaying and immediately flanking genomic regions.

## Change directory:
cd /home/ivan/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/

## Create folders:
mkdir -p outliers
mkdir -p outliers/data_outliers

## Loop over species:
for sp in Pcr Pfi Pga Pli Ppi Pva
do

## Loop over cutoffs:
#for c in 0.9900 0.9925 0.9950 0.9975 0.9990
for c in 0.9925
do

## Status:
echo "Now processing $sp"!

## Intersect BayPass SNPs with Tajima's D outlier windows:
bedtools intersect \
	-a baypass/0326/data_outliers/${sp}_median_outlier_SNPs_${c}.bed \
	-b TajimasD/0326/data_outliers/${sp}_outlier_windows_TajimasD_w5000_${c}.bed \
	-wa -u > outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_${c}.bed

## Expand SNPs 10 kb on both flanks:
bedtools slop \
	-i outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_${c}.bed \
	-g refgenomes/PodMur_1.0_GCF_004329235.1.chrom.sizes.chromnumbers \
	-b 10000 \
	> outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_10kb_${c}.bed

## Map to genes:
## Output columns will be: SNP_chr, SNP_start, SNP_end, gene_chr, gene_start, gene_end, gene, gene_ID
bedtools intersect \
	-a outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_10kb_${c}.bed \
	-b baypass/annotations/GCF_004329235.1_PodMur_1.0_genomic_patched_candidate_region_genes.bed -wa -wb \
	> outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_10kb_outlier_genes_${c}.bed

## Extract unique genes:
cut -f1,7,8 outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_10kb_outlier_genes_${c}.bed \
	| sort -u \
	> outliers/data_outliers/${sp}_SNPs_in_TajimasD_windows_10kb_outlier_genes_unique_${c}.txt

## Close loops:
done
done

## End of script.
