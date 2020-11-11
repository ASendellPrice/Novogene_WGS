#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --array=1-137:1
#SBATCH --job-name=filter_reads
#SBATCH --output=filter_reads.log
#SBATCH --error=filter_reads.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# CONDUCT FILTERING OF WGS RAW READS
# A. Sendell-Price, Sept 2020
#
# This script will apply the following filtering steps to raw
# sequencing reads received from Novogene:
# 1) Illumina adapters will be removed using trimmomatic
# 2) Duplicates will be removed using FastUniq

#########################################################################

# STEP 1:
# We will need to create a text file specifying the name of samples
# we want to process and the directory that raw reads are stored in,
# using 3 samples as an exaple, this will look like this:

# NOR140 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR141 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR142 /data/zool-zost/Novogene/NorfolkIsland/NOR142

Sample_List=/data/zool-zost/Novogene/sample.list.txt

# STEP 4:
# Use slurm array task ID to alocate sample name and directory
SAMPLE_NAME=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $2}')

# STEP 5:
#move into sample directory
cd $SAMPLE_DIRECTORY

# STEP 3:
#Load java
module load java/1.8.0

# STEP 4:
# Define path to trimmomatic, FastUniq and file specifying input for
# FastUniq:
trimmomatic_path=/data/zool-zost/BIN/Trimmomatic-0.39/trimmomatic-0.39.jar
fastuniq_path=/data/zool-zost/BIN/FastUniq/source/fastuniq
FastUniq_input=/data/zool-zost/Novogene/Scripts/fastqlist.txt
#Note: FastUniq_input is a simple text file that looks like this:
#AdapterTrimmed_forward_paired.fq
#AdapterTrimmed_reverse_paired.fq

# STEP 5:
# Set up for loop to conduct filtering for each sample
# 1..193 <- this will need to be amended depending on number of samples
# included in sample.list.txt
for ReadPair in `ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4 -d'_'`
do

    #Use trimmomatic to remove reads that dont meet quality threshold
    echo "Trimming reads ..."

    java -jar $trimmomatic_path PE \
    ${ReadPair}_1.fq.gz \
    ${ReadPair}_2.fq.gz \
    AdapterTrimmed_forward_paired.fq \
    AdapterTrimmed_forward_unpaired.fq \
    AdapterTrimmed_reverse_paired.fq \
    AdapterTrimmed_reverse_unpaired.fq \
    ILLUMINACLIP:/data/zool-zost/BIN/Trimmomatic-0.39/adapters/NovogeneIlluminaAdapters.fa:2:8:10 \
    HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:100

    #remove unpaired fastq files produced by Trimmomatic (these
    #wont be needed)
    rm *unpaired.fq

    #Run fastuniq
    echo "Removing diplicates reads ..."

    $fastuniq_path \
    -i $FastUniq_input \
    -t q \
    -o Dupfiltered_AdapterTrimmed_${ReadPair}_1.fq \
    -p Dupfiltered_AdapterTrimmed_${ReadPair}_2.fq

    #Compress dupfiltered and adapter trimmed fastq files
    gzip Dupfiltered_AdapterTrimmed_${ReadPair}*.fq -f

    #Remove intermediate fastq files
    rm *.fq

done
