#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-137:1
#SBATCH --time=10:00:00
#SBATCH --job-name=CreateSampleBams
#SBATCH --output=CreateSampleBams.log
#SBATCH --error=CreateSampleBams.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# Merge multiple bams per sample into a single sample bam file          #
# A. Sendell-Price, Sept 2020                                           #
#                                                                       #
# This script will take the bam files generated for each lane and merge #
# into a single bam file per sample.                                    #
#                                                                       #
# Note: requires output from script                                     #
# "/data/zool-zost/Novogene/Scripts/Map_sample_reads.sh"                #
# This script requires the number of jobs to match the number of lines  #
# in file sample.list.txt (update #SBATCH --array=1-137:1)              #                             
#########################################################################

#Load required modules
module load samtools

#Get sample information from sample list
SAMPLE_LIST=sample.list.txt
SAMPLE_NO=$SLURM_ARRAY_TASK_ID
SAMPLE_NAME=$(cat $SAMPLE_LIST | head -n $SAMPLE_NO | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $SAMPLE_LIST | head -n $SAMPLE_NO | tail -1 | awk {'print $2}')

#Move into sample directory
cd $SAMPLE_DIRECTORY

#Create list of sample bam files
ls *.bam > bam.list

#Merge bams using samtools merge
samtools merge -b bam.list ${SAMPLE_NAME}.bam

#Index bam file
samtools index ${SAMPLE_NAME}.bam

#remove bam list
rm bam.list
