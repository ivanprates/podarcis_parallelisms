module load bioinfo-tools BEDTools/2.29.2

#first collapse FST outliers in 20 kb windows
bedtools merge -i Outliers_IT_5kb2.5_FST.bed -c 1 -o count -d 20000 > Outliers_IT_5kb2.5_FST_merged_d10kb.bed
wc -l Outliers_IT_5kb2.5_FST_merged_d10kb.bed
bedtools merge -i Outliers_SA_5kb2.5_FST.bed -c 1 -o count -d 20000 > Outliers_SA_5kb2.5_FST_merged_d10kb.bed
wc -l Outliers_SA_5kb2.5_FST_merged_d10kb.bed

#outliers shared between FST IT and SA
bedtools intersect -a Outliers_IT_5kb2.5_FST_merged_d10kb.bed -b Outliers_SA_5kb2.5_FST_merged_d10kb.bed > Outliers_ITSA_5kb2.5_FST_merged_d10kb.bed
wc -l Outliers_ITSA_5kb2.5_FST_merged_d10kb.bed

###then cross-check with RAD-outliers
###first merge RAD outliers
bedtools merge -i Outliers_RAD_IT_new.bed -c 1 -o count -d 1 > Outliers_RAD_IT_new_merged.bed
wc -l Outliers_RAD_IT_new_merged.bed
bedtools merge -i Outliers_RAD_SA_new.bed -c 1 -o count -d 1 > Outliers_RAD_SA_new_merged.bed
wc -l Outliers_RAD_SA_new_merged.bed

bedtools intersect -wa -a Outliers_IT_5kb2.5_FST_merged_d10kb.bed -b Outliers_RAD_IT_new_merged.bed > Outliers_FSTRAD_IT_intersect_wa.bed
bedtools merge -i Outliers_FSTRAD_IT_intersect_wa.bed -d 1 > Outliers_FSTRAD_IT_intersect_wa_merged.bed
wc -l Outliers_FSTRAD_IT_intersect_wa_merged.bed
bedtools intersect -wa -a Outliers_SA_5kb2.5_FST_merged_d10kb.bed -b Outliers_RAD_SA_new_merged.bed > Outliers_FSTRAD_SA_intersect_wa.bed
bedtools merge -i Outliers_FSTRAD_SA_intersect_wa.bed -d 1 > Outliers_FSTRAD_SA_intersect_wa_merged.bed
wc -l Outliers_FSTRAD_SA_intersect_wa_merged.bed

###last, cross-check with Tajima's D outliers
###first merge Tajima outliers
bedtools merge -i Outliers_TajimasD_IT.bed -c 1 -o count -d 1 > Outliers_TajimasD_IT_merged.bed
wc -l Outliers_TajimasD_IT_merged.bed
bedtools merge -i Outliers_TajimasD_SA.bed -c 1 -o count -d 1 > Outliers_TajimasD_SA_merged.bed
wc -l Outliers_TajimasD_SA_merged.bed

bedtools intersect -wa -a Outliers_FSTRAD_IT_intersect_wa_merged.bed -b Outliers_TajimasD_IT_merged.bed > Outliers_FSTRADTajima_IT_intersect_wa.bed
bedtools merge -i Outliers_FSTRADTajima_IT_intersect_wa.bed -d 1 > Outliers_FSTRADTajima_IT_intersect_wa_merged.bed
wc -l Outliers_FSTRADTajima_IT_intersect_wa_merged.bed
bedtools intersect -wa -a Outliers_FSTRAD_SA_intersect_wa_merged.bed -b Outliers_TajimasD_SA_merged.bed > Outliers_FSTRADTajima_SA_intersect_wa.bed
bedtools merge -i Outliers_FSTRADTajima_SA_intersect_wa.bed -d 1 > Outliers_FSTRADTajima_SA_intersect_wa_merged.bed
wc -l Outliers_FSTRADTajima_SA_intersect_wa_merged.bed

###finally, cross-check three-way outliers between IT and SA
bedtools intersect -a Outliers_FSTRADTajima_IT_intersect_wa_merged.bed -b Outliers_FSTRADTajima_SA_intersect_wa_merged.bed > Outliers_FSTRADTajima_ITSA_intersect_wa_merged.bed
more Outliers_FSTRADTajima_ITSA_intersect_wa_merged.bed 
