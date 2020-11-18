#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=52:00:00
#SBATCH --array=1-10:1
#SBATCH --job-name=Batch1_Pipeline
#SBATCH --output=Batch1.%A_%a.out
#SBATCH --error=Batch1.%A_%a.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#######################################################################################################
# SET SAMPLE NAME AND DIRECTORY
# This pipeline requires the number of jobs to match the number of lines in file batch1.txt (see below) 
# (update #SBATCH --array=1-137:1)
#######################################################################################################

# STEP 1:
# We will need to create a text file specifying the name of samples we want to process and the directory
# that raw reads are stored in, using 3 samples as an exaple, this will look like this:
# NOR140 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR141 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR142 /data/zool-zost/Novogene/NorfolkIsland/NOR142

#Set directory to that list
Sample_List=/data/zool-zost/Novogene/batch1.txt

# STEP 2:
# Use slurm array task ID to alocate sample name and directory
SAMPLE_NAME=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $2}')

# STEP 3:
#move into sample directory
cd $SAMPLE_DIRECTORY


#######################################################################################################
# CONDUCT FILTERING OF WGS RAW READS
# We will use FastP - a tool designed to provide fast all-in-one preprocessing for FastQ files
# see https://github.com/OpenGene/fastp for more info.
#######################################################################################################

# STEP 1:
# Define path to Fastp:
FASTP=/data/zool-zost/BIN/fastp

# STEP 2:
# Set up for loop to conduct filtering for each read pair
for ReadPair in `ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4 -d'_'`
do
  
  #Use Fastp to conduct automated filtering of fastq files
  #Note: based on initial test we will trim the first 10bp from start of each read
  $FASTP \
  -i ${ReadPair}_1.fq.gz \
  -o Filtered_${ReadPair}_1.fq.gz \
  -I ${ReadPair}_2.fq.gz \
  -O Filtered_${ReadPair}_2.fq.gz \
  --trim_front1 10 \
  --trim_front2 10

  #Remove .json file as not needed
  rm fastp.json

  #Rename QC report and move to fastp qc report folder
  mv fastp.html /data/zool-zost/Novogene/fastp_QC_reports/${ReadPair}.html

done


#######################################################################################################
# MAP FILTERED READS TO REFERENCE ASSEMBLY
#######################################################################################################

# STEP 1:
# Load required modules in ARC
module load java/1.8.0
module load samtools
module load bowtie2

# STEP 2:
# Specify "plate" name for sequencing run
PLATE=Novogene_Batch1_WGS_2020

# STEP 3:
# Specify path to reference genome (cornetti assembly) including file prefix "ZOLAv0"
GENOME_DB=/data/zool-zost/Ref_Genome/ZOLAv0

# STEP 4:
# For each pair of reads conduct mapping using bowtie2
for ReadPair in `ls Filtered_${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5 -d'_'`
do
  	bowtie2 -x $GENOME_DB \
    -1 ${ReadPair}_1.fq.gz \
    -2 ${ReadPair}_2.fq.gz \
    --rg-id $ReadPair --rg SM:$SAMPLE_NAME --rg LB:$PLATE \
    --rg PU:$PLATE --rg PL:illumina 2> ${ReadPair}.log  | \
    samtools view -bhS - | \
    samtools sort - > ${ReadPair}.bam
done


#######################################################################################################
# MERGE SAMPLE BAMS INTO A SINGLE FILE
# Currently we have a bam file for each read pair per sample. We will now merge these to create a
# single bam file per sample
#######################################################################################################

# Count number of ma files per sample
BAM_Count=$(ls Filtered_${SAMPLE_NAME}*.bam | wc -l)

# If number of bams is greater than 1 then merge bams into a single file using samtools merge
# else, rename single bam file
if [ $BAM_Count -gt 1 ]
then
  samtools merge ${SAMPLE_NAME}.bam Filtered_${SAMPLE_NAME}*.bam 
  rm Filtered_${SAMPLE_NAME}*.bam
else
  mv Filtered_${SAMPLE_NAME}*.bam ${SAMPLE_NAME}.bam
fi

# Sort bam file with samtools sort - required for genotyping with ANGSD
samtools sort ${SAMPLE_NAME}.bam > ${SAMPLE_NAME}.sorted.bam
rm ${SAMPLE_NAME}.bam

#######################################################################################################
# CONDUCT SAMPLE GENOTYPING
# Use GATK haployype calling to generate an intermediate sample gVCF file containing called genotypes 
# at both variable and non-variable sites. We will run haplotype calling using the following settings:
# EMIT_ALL_CONFIDENT_SITES = output both variant and non-variant sites
# BP_RESOLUTION = output invariant sites bp by bp rather than as blocks.
#######################################################################################################

# STEP 1: 
# Set path to GATK Jar file
GATK_JAR=/data/zool-zost/BIN/gatk-4.1.9.0/gatk

#STEP 2:
# Set path to reference genome assembly (silvereye v.1)
# This is the Cornetti assembly
FASTA=/data/zool-zost/Ref_Genome/GCA_001281735.1_ASM128173v1_genomic.fna

#STEP 3:
# Call geotypes with GATK haplotype caller.
$GATK_JAR HaplotypeCaller \
-R $FASTA \
-I ${SAMPLE_NAME}.sorted.bam \
--emit-ref-confidence BP_RESOLUTION \
--output-mode EMIT_ALL_CONFIDENT_SITES \
-O ${SAMPLE_NAME}.g.vcf.gz

