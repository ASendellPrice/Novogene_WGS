#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --array=1-137:1
#SBATCH --job-name=Map_Reads
#SBATCH --output=Map_Reads.log
#SBATCH --error=Map_Reads.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# Mapping filtered reads to reference assembly                          #
# A. Sendell-Price, Sept 2020                                           #
#                                                                       #
# This script will map filtered Novogene WGS sequencing reads to the    #
# cornetti et al. silvereye reference genome version 1.                 #
#                                                                       #
# Note: requires output from script                                     #
# "/data/zool-zost/Novogene/Scripts/FilterReads.sh"                     #
# This script requires the number of jobs to match the number of lines  #
# in file sample.list.txt (update #SBATCH --array=1-137:1)              #
#########################################################################

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
# Specify path to sample list - this will be the same file used in
# FilterReads.sh
Sample_List=/data/zool-zost/Novogene/sample.list.txt

# STEP 5:
#Use slurm array task ID to alocate sample name and directory
SAMPLE_NAME=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $2}')

# STEP 6:
#move into sample directory
cd $SAMPLE_DIRECTORY

# STEP 7:
#For each pair of reads conduct mapping using bowtie2
for ReadPair in `ls Dupfiltered_AdapterTrimmed_${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5,6 -d'_'`
do
  	bowtie2 -x $GENOME_DB \
    -1 ${ReadPair}_1.fq.gz \
    -2 ${ReadPair}_2.fq.gz \
    --rg-id $ReadPair --rg SM:$SAMPLE_NAME --rg LB:$PLATE \
    --rg PU:$PLATE --rg PL:illumina 2> ${ReadPair}.log  | \
    samtools view -bhS - | \
    samtools sort - > ${ReadPair}.bam
done
