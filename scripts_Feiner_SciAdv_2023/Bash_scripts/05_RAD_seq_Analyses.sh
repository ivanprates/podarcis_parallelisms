#after demultiplexing, RAD-seq data was processed according to the stacks pipeline

module load bioinfo-tools
module load samtools/1.12
module load bwa/0.7.17

#first reads were aligned to the reference genome
#this step was done for all samples; here shown for an example
bwa mem -t 12 /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
      CA51.1.fq.gz CA51.2.fq.gz | samtools sort --threads 12 -O bam -o CA51.bam

#then stacks with population was used with popmap_900ind, which is a list of all samples and their population
module load Stacks/2.59

gstacks -I ./all_bams/ -O ./stacksFiles/ -M popmap_900ind -t 40

populations -P ./stacksFiles/ -M popmap_900ind -t 40 -R 0.95 \
--min-maf 0.03 --vcf --plink

#files were filtered for genotype rate, missing data and LD, keeping only chromosomes:
cd stacksFiles
plink -file populations.plink -aec --set-missing-var-ids @:# --geno 0.5 --mind 0.5 \
--chr CM014743.1 CM014744.1 CM014745.1 CM014746.1 CM014747.1 CM014748.1 CM014749.1 CM014750.1 CM014751.1 CM014752.1 CM014753.1 CM014754.1 CM014755.1 CM014756.1 CM014757.1 CM014758.1 CM014759.1 CM014760.1 CM014761.1 \
-indep-pairwise 1 kb 1 0.8 -out ldkb1r0.8
plink -file populations.plink -aec --set-missing-var-ids @:# --geno 0.5 --mind 0.5 \
--chr CM014743.1 CM014744.1 CM014745.1 CM014746.1 CM014747.1 CM014748.1 CM014749.1 CM014750.1 CM014751.1 CM014752.1 CM014753.1 CM014754.1 CM014755.1 CM014756.1 CM014757.1 CM014758.1 CM014759.1 CM014760.1 CM014761.1 \
-exclude ldkb1r0.8.prune.out -out populations.plink_pruned -make-bed

#check for missing genotypes; this creates .imissi file that contains missing genotype rate per sample
plink -bfile populations.plink_pruned --missing --aec --set-missing-var-ids @:# --geno 0.5 --mind 0.5

#PCA
plink -bfile populations.plink_pruned --pca --out populations_PCA --aec --set-missing-var-ids @:# --geno 0.5 --mind 0.5

#run structure analysis
python $fastStructure_ROOT/structure.py -K 2 --cv=10 --input=populations.plink_pruned --output=structure_result

#check how many SNPs are in candidate region (cand_region.bed)
plink -bfile populations.plink_pruned --aec --set-missing-var-ids @:# --geno 0.5 --mind 0.5 --recode vcf
vcftools --vcf plink.vcf --bed cand_region.bed --out plink_cand --recode --recode-INFO-all
paste plink_cand.recode.vcf | grep -v "#" | awk '{ print $1, $2}' > SNPs_candreg

###after this the R_script RADseq_SplitDatasetITSA.R should be run. This produces lineage-specific (IT, SA, hybrid) data files and co-variate files (with sex and PCs of genetic variations)

###Thereafter, GWAS analyses are run using plink for the IT lineage by taking sex and PC1 as co-variate:
plink -bfile out --aec --allow-no-sex --set-missing-var-ids @:# --linear --adjust --covar MyCov --covar-number 1-2

#and for the SA lineage by taking sex as co-variate (PC1 is confounded with greenness)
plink -bfile out --aec --allow-no-sex --set-missing-var-ids @:# --linear --adjust --covar MyCov --covar-number 1

###after this the R_script RADseq_Export_Outliers.R is run to export outliers.

