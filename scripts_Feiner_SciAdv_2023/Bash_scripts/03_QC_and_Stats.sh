#This script gathers a number of statistics on LD decay, missing genotyping data and read-coverage and SNP-densities in the candidate region.

module load bioinfo-tools 
module load BEDTools/2.29.2 

#obtaining the number of SNPs per 10kb-window before and after filtering:
bedtools makewindows -g /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.genomeFile.txt -w 10000 > windows.bed
bedtools coverage -a windows.bed -b ../Filtering/GATK_all_SNP_final_cand20mb.recode.vcf -counts > GATK_all_SNP_filtered_pass_cand_plusminus20mb_coverage_10kb
bedtools coverage -sorted -a windows.bed -b ../Filtering/GATK_all_SNP_final_prefiltered_cand20mb.vcf.recode.vcf -counts > GATK_all_SNP_filtered_cand_plusminus20mb_coverage_10kb

#obtaining the raw read depth (DP; pre-filtering) per sample per SNP:
grep -v '#' ../Filtering/GATK_all_SNP_final_prefiltered_cand.vcf.recode.vcf | cut -f 10-159 | awk '{for(i=1;i<=NF;++i) cols[i]} {for(i=1;i<=NF;i++) if(i in cols) {split($i,t,":"); $i=t[3]}}1' > tmp1
grep -v '##' ../Filtering/GATK_all_SNP_final_prefiltered_cand.vcf.recode.vcf | head -1 | cut -f 10-159 > IDs
cat IDs tmp1 > tmp2
grep -v '##' ../Filtering/GATK_all_SNP_final_prefiltered_cand.vcf.recode.vcf | cut -f 2 > pos
paste pos tmp2 > Read_Depth_Per_Ind
sed -i 's/ /\t/g' Read_Depth_Per_Ind
rm tmp* IDs pos

# the resulting files were imported into R and plotted with the R-script Coverage_stats.R


###TE density in and around candidate region
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_rm.out.gz
gunzip GCF_004329235.1_PodMur_1.0_rm.out.gz

awk '{if ($5 == "NC_041323.1") print $0;}' GCF_004329235.1_PodMur_1.0_rm.out | awk '{print $6,$7,$11}' > Chr12_TEs

mkdir Chr12_10kb
cp Chr12_100kb/Chr12_TEs Chr12_10kb/.
cd Chr12_10kb

for (( i=1; i <= 6110; i++ ))
do
grep " $i[0-9][0-9][0-9][0-9] " Chr12_TEs > $(($i+1))0000
done

wc -l *000 > Chr12_binning_summary_10kb

# the resulting files were imported into R and plotted with the R-script TE_density.R


###LD-decay for each chromosome

module load bioinfo-tools
module load plink/1.90b4.9
module load python/2.7.11

#run plink for every chromosome (here chr. 1 is shown as example) and summarize decay for each chromosome using python script 'ld_decay_calc.py' (downloaded from github: https://github.com/speciationgenomics/scripts/blob/master/ld_decay_calc.py)
plink --vcf ../Filtering/GATK_all_SNP_final.vcf -aec --set-missing-var-ids @:# --chr CM014743.1 --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 50 --ld-window-r2 0 --out LD_decay_chr1 --mind 0.5 --geno 0.5
python ld_decay_calc.py -i LD_decay_chr1.ld.gz -o LD_decay_chr1

# the resulting files were imported into R and plotted with the R-script R_LD_decay_final.R


###pairwise LD in and around candidate region

module load bioinfo-tools
module load plink/1.90b4.9

plink --vcf ../Filtering/GATK_all_SNP_final_cand21_22.recode.vcf --allow-extra-chr --set-missing-var-ids @:# --r2 --ld-window 100000 --ld-window-kb 1000 --ld-window-r2 0 --out 1mb_pairwise_ld

# the resulting files were imported into R and plotted with the R-script LD_heatmap.R
