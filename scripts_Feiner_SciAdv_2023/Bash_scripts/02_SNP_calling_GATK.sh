#This script uses GATK to call SNPs. 
#To speed up the process, mpileup was run for each chromosome separately

module load bioinfo-tools
module load bwa/0.7.17
module load picard/2.23.4
module load GATK/4.2.0.0
module load FastQC/0.11.9
module load trimmomatic/0.36
module load samtools/1.14
module load plink/1.90b4.9

mkdir /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24

#Run HaplotypeCaller
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014743.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014743.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014744.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014744.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014745.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014745.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014746.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014746.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014747.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014747.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014748.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014748.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014749.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014749.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014750.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014750.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014751.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014751.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014752.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014752.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014753.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014753.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014754.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014754.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014755.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014755.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014756.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014756.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014757.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014757.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014758.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014758.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014759.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014759.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014760.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014760.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L CM014761.1 -ERC GVCF -O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_CM014761.g.vcf.gz &
gatk --java-options "-Xmx6g" HaplotypeCaller -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -I mapping/AS24_marked_renamed.bam \
-L /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic_nonChr.list -ERC GVCF \
-O /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24/AS24_Rest.g.vcf.gz &
wait

cd /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Mappings/HaplotypeCaller_Output/AS24

#index final g.vcf file:
gatk IndexFeatureFile -I AS24_CM014743.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014744.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014745.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014746.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014747.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014748.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014749.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014750.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014751.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014752.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014753.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014754.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014755.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014756.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014757.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014758.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014759.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014760.g.vcf.gz
gatk IndexFeatureFile -I AS24_CM014760.g.vcf.gz
gatk IndexFeatureFile -I AS24_Rest.g.vcf.gz

#then for each chromosome run CombineGVCFs, GenotypeGVCFs and VariantFiltration

mkdir Tmp_dir_comb_CM014743
mkdir Tmp_dir_comb_CM014743_2

gatk --java-options "-Xmx125g" CombineGVCFs -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
-V input_CM014743.list -O CM014743.g.vcf.gz --tmp-dir Tmp_dir_comb_CM014743

gatk --java-options "-Xmx125g" GenotypeGVCFs -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
-V CM014743.g.vcf.gz -O CM014743.vcf.gz --tmp-dir Tmp_dir_comb_CM014743_2

rm -rf Tmp_dir_comb_CM014743
rm -rf Tmp_dir_comb_CM014743_2

gatk --java-options "-Xmx125g" SelectVariants -V CM014743.vcf.gz -select-type SNP -O CM014743_SNP.vcf.gz
gatk IndexFeatureFile -I CM014743_SNP.vcf.gz
gatk --java-options "-Xmx125g" SelectVariants -V CM014743.vcf.gz -select-type INDEL -O CM014743_INDEL.vcf.gz
gatk IndexFeatureFile -I CM014743_INDEL.vcf.gz

gatk --java-options "-Xmx125g" VariantFiltration -R /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna -V CM014743_SNP.vcf.gz \
--mask CM014743_INDEL.vcf.gz \
--mask-extension 10 \
--mask-name InDel \
--filter-expression "FS > 60.0" --filter-name "FisherStrand" \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "MQ < 40.0" --filter-name "MapQual40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--filter-expression "DP < 1192"  --filter-name "minDepth" \
--filter-expression "DP > 4768"  --filter-name "maxDepth" \
-O CM014743_SNP_filtered.vcf.gz

#Lastly, combine all chromosome files into one
java -jar $PICARD_HOME/picard.jar GatherVcfs \
I=CM014743_SNP_filtered.vcf.gz \
I=CM014744_SNP_filtered.vcf.gz \
I=CM014745_SNP_filtered.vcf.gz \
I=CM014746_SNP_filtered.vcf.gz \
I=CM014747_SNP_filtered.vcf.gz \
I=CM014748_SNP_filtered.vcf.gz \
I=CM014749_SNP_filtered.vcf.gz \
I=CM014750_SNP_filtered.vcf.gz \
I=CM014751_SNP_filtered.vcf.gz \
I=CM014752_SNP_filtered.vcf.gz \
I=CM014753_SNP_filtered.vcf.gz \
I=CM014754_SNP_filtered.vcf.gz \
I=CM014755_SNP_filtered.vcf.gz \
I=CM014756_SNP_filtered.vcf.gz \
I=CM014757_SNP_filtered.vcf.gz \
I=CM014758_SNP_filtered.vcf.gz \
I=CM014759_SNP_filtered.vcf.gz \
I=CM014760_SNP_filtered.vcf.gz \
I=CM014761_SNP_filtered.vcf.gz \
O=GATK_all_SNP_filtered.vcf.gz

cp GATK_all_SNP_filtered.vcf.gz /proj/snic2020-16-107/nobackup/private/Analysis_149samples/GATK/.
cd /proj/snic2020-16-107/nobackup/private/Analysis_149samples/GATK
gatk IndexFeatureFile -I GATK_all_SNP_filtered.vcf.gz
bcftools index GATK_all_SNP_filtered.vcf.gz

vcftools --gzvcf GATK_all_SNP_filtered.vcf.gz --bed cand_region.bed --out GATK_all_SNP_filtered_cand.vcf --recode --recode-INFO-all #do this do get coverage-stats for un-filtered read data

bcftools view GATK_all_SNP_filtered.vcf.gz -f 'PASS,.' -Ob -o GATK_all_SNP_filtered_pass.bcf
bcftools index GATK_all_SNP_filtered_pass.bcf

#further filtering:
#for bi-allelic SNPs:
bcftools view -m2 -M2 -v snps ../GATK_all_SNP_filtered_pass.bcf -Ob -o GATK_all_SNP_filtered_pass_biallelic.bcf
bcftools index GATK_all_SNP_filtered_pass_biallelic.bcf

#for minor allele frequency:
bcftools view -q 0.03:minor GATK_all_SNP_filtered_pass_biallelic.bcf -Ob -o GATK_all_SNP_filtered_pass_biallelic_maf0.03.bcf
bcftools index GATK_all_SNP_filtered_pass_biallelic_maf0.03.bcf

#for genotype quality (GQ) below 20:
bcftools filter -S . -e 'FMT/GQ<20' GATK_all_SNP_filtered_pass_biallelic_maf0.03.bcf -Ob -o GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf
bcftools index GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf

#delete SNPs within 10 bp of each other:
bcftools +prune -w 10bp -n 1 -N rand GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf -Ob -o GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf
bcftools index GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf

#delete SNPs with more than 50% missing genotype data
bcftools filter -e 'F_MISSING > 0.5' GATK_all_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf -Ov -o GATK_all_SNP_final.vcf
bcftools stats GATK_all_SNP_final.vcf > GATK_all_SNP_final.vchk #this returns stats on 20,928,004 SNPs)

#extract candidate region:
vcftools --vcf GATK_all_SNP_final.vcf --bed ../cand_region.bed --out GATK_all_SNP_final_cand --recode --recode-INFO-all
vcftools --vcf GATK_all_SNP_final.vcf --bed ../cand_region_plusminus1mb.bed --out GATK_all_SNP_final_cand1mb --recode --recode-INFO-all
vcftools --vcf GATK_all_SNP_final.vcf --bed ../cand_region_plusminus20mb.bed --out GATK_all_SNP_final_cand20mb --recode --recode-INFO-all
vcftools --vcf GATK_all_SNP_final.vcf --bed ../cand_region_21_22.bed --out GATK_all_SNP_final_cand21_22 --recode --recode-INFO-all

#run PCA analysis
plink --vcf GATK_all_SNP_final.vcf --pca --out populations_PCA --aec
# the resulting files were imported into R and plotted with the R-script Plot_PCA_WGS.R

#same for pre-filtered SNPs
vcftools --gzvcf ../GATK_all_SNP_filtered.vcf.gz --bed ../cand_region.bed --out GATK_all_SNP_final_prefiltered_cand.vcf --recode --recode-INFO-all
vcftools --gzvcf ../GATK_all_SNP_filtered.vcf.gz --bed ../cand_region_plusminus1mb.bed --out GATK_all_SNP_final_prefiltered_cand1mb.vcf --recode --recode-INFO-all
vcftools --gzvcf ../GATK_all_SNP_filtered.vcf.gz --bed ../cand_region_plusminus20mb.bed --out GATK_all_SNP_final_prefiltered_cand20mb.vcf --recode --recode-INFO-all

#note: the file GATK_all_SNP_final.vcf will be used as input for FST-outlier scans and Tajima's D analyses (see script 04)
