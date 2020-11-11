#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --array=1-137:1
#SBATCH --job-name=InitalQC
#SBATCH --output=InitalQC.log
#SBATCH --error=InitalQC.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# CONDUCT INITIAL QUALITY CHECK OF SEQUENCING READS
# A. Sendell-Price, Sept 2020
#
# This script will perform initial quality checking of raw sequencing reads
# received from Novogene using fastQC
#
# This script requires the number of jobs to match the number of lines
# in file sample.list.txt (update #SBATCH --array=1-137:1)
#########################################################################

# STEP 1:
# Define path to fastQC:
FASTQC=/data/zool-zost/BIN/FastQC_2GigHack/fastqc

# STEP 2:
# Load java module needed to run FastQC
module load java/1.8.0

# STEP 3:
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

for FastQ_file in `ls ${SAMPLE_NAME}_*.fq.gz`
do

  #Pass sequencing reads to FastQC
  zcat ${FastQ_file} | $FASTQC stdin
  mv stdin_fastqc.html ${FastQ_file}.html
  rm stdin_fastqc.zip
  
done
