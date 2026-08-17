#this script was run for each of the 149 samples;
#each sample had it's own folder including the raw fastq files
#here is the example of sample AS24

module load bioinfo-tools
module load bwa/0.7.17
module load picard/2.23.4
module load GATK/4.2.0.0
module load FastQC/0.11.9
module load trimmomatic/0.39
module load samtools/1.14

# new directories
mkdir fastqcRaw
mkdir fastqcTrim
mkdir trimReads
mkdir mapping

# QC
fastqc -o fastqcRaw -f fastq -t 16 AS24_1.fq.gz AS24_2.fq.gz

# trim
java -jar $TRIMMOMATIC_ROOT/trimmomatic-0.39.jar PE -threads 16 -phred33 \
AS24_1.fq.gz AS24_2.fq.gz \
trimReads/AS24_trim_1.fq.gz trimReads/AS24_trim_unpair1.fq.gz \
trimReads/AS24_trim_2.fq.gz trimReads/AS24_trim_unpair2.fq.gz \
ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70

# QC after trim
fastqc -o fastqcTrim/ -f fastq -t 16 trimReads/*

# fastqc results before and after trimming were checked manually one by one by opening the results files in browser. E.g.:
xdg-open fastqcRaw/LE24_1_fastqc.html

# mapping to the reference genome; here against the Spanish reference genome, but this was also repeated for the newly assembled genomes (SMA and CDS)
bwa mem -M -t 16 /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
trimReads/AS24_trim_1.fq.gz trimReads/AS24_trim_2.fq.gz > mapping/AS24_pair.sam

bwa mem -M -t 16 /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
trimReads/AS24_trim_unpair1.fq.gz > mapping/AS24_unpair_1.sam

bwa mem -M -t 16 /crex/proj/snic2020-6-73/Projects/Wallie_Genomics/Genomes/Pmu_Spanish/GCA_004329235.1_PodMur_1.0_genomic.fna \
trimReads/AS24_trim_unpair2.fq.gz > mapping/AS24_unpair_2.sam

# sam to bam
java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups INPUT=mapping/AS24_pair.sam OUTPUT=mapping/AS24_pair.bam SORT_ORDER=coordinate \
RGLB=1708 RGPL=ILLUMINA RGPU=386 RGSM=AS24

java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups INPUT=mapping/AS24_unpair_1.sam OUTPUT=mapping/AS24_unpair_1.bam SORT_ORDER=coordinate \
RGLB=1708 RGPL=ILLUMINA RGPU=386 RGSM=AS24

java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups INPUT=mapping/AS24_unpair_2.sam OUTPUT=mapping/AS24_unpair_2.bam SORT_ORDER=coordinate \
RGLB=1708 RGPL=ILLUMINA RGPU=386 RGSM=AS24

# merge bam
java -jar $PICARD_HOME/picard.jar MergeSamFiles INPUT=mapping/AS24_pair.bam INPUT=mapping/AS24_unpair_1.bam INPUT=mapping/AS24_unpair_2.bam \
OUTPUT=mapping/AS24.bam SORT_ORDER=coordinate ASSUME_SORTED=true

# index bam
java -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=mapping/AS24.bam

# remove intermediate files
rm mapping/AS24_*

#check the mapping stats:
samtools flagstat mapping/AS24.bam > mapping/AS24_MappingStats

#mark duplicates
java -jar $PICARD_HOME/picard.jar MarkDuplicates INPUT=mapping/AS24.bam OUTPUT=mapping/AS24_marked.bam METRICS_FILE=mapping/AS24_marked.metrics

#Index .bam file again
java -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=mapping/AS24_marked.bam

# remove intermediate files
rm mapping/AS24.ba*

#add meaningful IDs
java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups INPUT=mapping/AS24_marked.bam OUTPUT=mapping/AS24_marked_renamed.bam RGID=AS24 RGLB=AS24-lib \
RGPL=ILLUMINA RGPU=AS24-01 RGSM=AS24

#index again
java -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=mapping/AS24_marked_renamed.bam

# remove intermediate files
rm mapping/AS24_marked.ba*
rm mapping/AS24.ba*
