module load bioinfo-tools
module load star/2.7.9a

#Raw reads were trimmed using trimmomatic; here one sample is shown as example
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE Sample1_R1.fastq.gz Sample1_R1.fastq.gz Sample1_F  Sample1_F_unpaired   Sample1_R  Sample1_R_unpaired  \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

#The reference genome was prepped
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir Reference_Genome --genomeFastaFiles Reference_Genome/GCF_004329235.1_PodMur_1.0_genomic.fna \
--sjdbGTFfile Reference_Genome/GCF_004329235.1_PodMur_1.0_genomic.gtf --sjdbGTFtagExonParentTranscript Parent

#Trimmed reads were mapped against the reference genome; here one sample is shown as example; only paired reads in which each read survived the trimming were used
STAR --runThreadN 16 --genomeDir Reference_Genome --readFilesIn Sample1_F Sample1_R --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts  --outFileNamePrefix STAR/Sample1

#summarize the individual count files into a counttable that contains counts from all samples:
cd STAR
paste *ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp
sed -e "1igene_name\t$(ls *ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-46 > ../raw_counts_matrix.txt #46 is because of 45 samples + 1 column for gene names
rm tmp

#make a file that contains positional data for each gene 
cd Reference_Genome
cat GCF_004329235.1_PodMur_1.0_genomic.gtf | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="gene"){split($9, a, "\""); print $1, $4, $5, a[2]}}' > Gene_positions


###Repeat the previous four steps with a custom gtf file that contains manually-curated gene models of the candidate region
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir Reference_Genome --genomeFastaFiles Reference_Genome/GCF_004329235.1_PodMur_1.0_genomic.fna \
--sjdbGTFfile Reference_Genome/Candidate_region_manually_annotated.gtf --sjdbGTFtagExonParentTranscript Parent
STAR --runThreadN 16 --genomeDir Reference_Genome --readFilesIn Sample1_F Sample1_R --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts  --outFileNamePrefix STAR_Cand/Sample1
cd STAR_Cand
paste *ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp
sed -e "1igene_name\t$(ls *ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-46 > ../raw_counts_matrix_cand.txt #46 is because of 45 samples + 1 column for gene names
rm tmp
cd Reference Genome
cat Candidate_region_manually_annotated.gtf | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="gene"){split($9, a, "\""); print $1, $4, $5, a[2]}}' > Gene_positions_cand

#the raw_counts_matrix* files were opened in excel and the triplicates per sample were manually collapsed into a single (average) measure (files 'raw_counts_matrix_RepCollapsed_rounded.csv' and 'raw_counts_matrix_Cand_RepCollapsed_rounded.csv')
###after this the R_script RNA_seq_analysis_QC_PCA.R was used to do QC and plot PCAs and Heatmaps;
###R_scripts DESeq2_analysis_BrownGreen_RepColl_Cand.R and DESeq2_analysis_BrownBlack_RepColl_Cand.R were used to run the DESeq analyses
