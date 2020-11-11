#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=52:00:00
#SBATCH --array=1-136:1
#SBATCH --job-name=genotyping-array
#SBATCH --output=genotyping-array.log
#SBATCH --error=genotyping-array.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# Conduct sample genotyping using GATK Haplotype caller                 #
# A. Sendell-Price, Sept 2020                                           #
#                                                                       #
# This script will generate an intermediate sample gVCF file containing #
# called genotypes at both variable and non-variable sites              #
#                                                                       #
# This script requires the number of jobs to match the number of lines  #
# in file sample.list.txt (update #SBATCH --array=1-137:1)              #
#########################################################################

# Load java module (needed to launch GATK)
module load java/1.8.0

# Store paths to necessary files/programs as variables
# The reference genome assembly (silvereye v.1)
FASTA=/data/zool-zost/Ref_Genome/GCA_001281735.1_ASM128173v1_genomic.fna

# Get sample name using slurm task ID
SAMPLE_NAME=$(cat Sample.list | head -n $SLURM_ARRAY_TASK_ID | tail -1 )

# Get path to sample bam file
BAM=/data/zool-zost/Novogene/

_Sample_Bams/${SAMPLE_NAME}.bam

# GATK Jar file
GATK_JAR=/data/zool-zost/BIN/gatk-4.1.9.0/gatk   

# Call geotypes with GATK haplotype caller.
$GATK_JAR --java-options "-Xmx4g" HaplotypeCaller \
-R $FASTA \
-I $BAM \
--emit-ref-confidence BP_RESOLUTION \
--output-mode EMIT_ALL_CONFIDENT_SITES \
-O ${SAMPLE_NAME}.g.vcf.gz
