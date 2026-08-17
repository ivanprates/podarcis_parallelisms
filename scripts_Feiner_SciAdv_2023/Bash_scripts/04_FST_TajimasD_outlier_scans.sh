#This script does FST-outlier scans and looks at patterns of genetic diversity using Tajima's D

module load bioinfo-tools
module load vcftools/0.1.16

###FST
vcftools --vcf GATK_all_SNP_final.vcf --weir-fst-pop Green_indiv_list_IT --weir-fst-pop Brown_indiv_list_IT --fst-window-size 5000 --fst-window-step 2500 --out FSTscan_GATK_20vs20_5kb2.5_IT
vcftools --vcf GATK_all_SNP_final.vcf --weir-fst-pop Green_indiv_list_SA --weir-fst-pop Brown_indiv_list_SA --fst-window-size 5000 --fst-window-step 2500 --out FSTscan_GATK_10vs10_5kb2.5_SA

#doing the same for single-SNPs for both lineages together (i.e., 30 vs 30)
vcftools --vcf GATK_all_SNP_final_new.vcf --weir-fst-pop Green_indiv_list_all --weir-fst-pop Brown_indiv_list_all --out FSTscan_GATK_30vs30
#using a custom script, these were then summarized for 100-SNP-windows (median values).

# the resulting files were imported into R and plotted with the R-script FSTscan_GATK_100SNPs_30vs30_Final_new_PlotCandRegion_FST_RAD_Tajima.R; export outlier to overlap-analysis with the R-script FSTscan_GATK_export_outliers_ITSA.R

###Tajima's D
vcftools --vcf GATK_all_SNP_final.vcf --keep Green_indiv_list_SA --TajimaD 5000 --out WG_5kb_TajimaD_green_SA
vcftools --vcf GATK_all_SNP_final.vcf --keep Brown_indiv_list_SA --TajimaD 5000 --out WG_5kb_TajimaD_brown_SA
vcftools --vcf GATK_all_SNP_final.vcf --keep Green_indiv_list_IT --TajimaD 5000 --out WG_5kb_TajimaD_green_IT
vcftools --vcf GATK_all_SNP_final.vcf --keep Brown_indiv_list_IT --TajimaD 5000 --out WG_5kb_TajimaD_brown_IT

# the resulting files were imported into R and plotted with the R-script FSTscan_GATK_100SNPs_30vs30_Final_new_PlotCandRegion_FST_RAD_Tajima.R; export outlier to overlap-analysis with the R-script Tajima_WG_ITSA_Final.R
