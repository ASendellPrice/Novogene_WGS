#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=52:00:00
#SBATCH --job-name=Merge_gVCFs
#SBATCH --output=Merge_gVCFs.log
#SBATCH --error=Merge_gVCFs.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#########################################################################
# Merge individual sample gVCFs into a single multi-sample gVCF         #
# A. Sendell-Price, Sept 2020                                           #
#                                                                       #
# Note: requires output from script                                     #
# "/data/zool-zost/Novogene/Scripts/Sample_Genotyping.sh".              #
# You will need to specify the sample gVCFs to be merged.               #
#########################################################################

# Load java module (needed to launch GATK)
module load java/1.8.0

# Store paths to necessary files/programs as variables
# The reference genome assembly (silvereye v.1)
FASTA=/data/zool-zost/Ref_Genome/GCA_001281735.1_ASM128173v1_genomic.fna

# GATK Jar file
GATK_JAR=/data/zool-zost/BIN/gatk-4.1.9.0/gatk  

# Combine sample gVCFs into single multi-sample gVCF
# Using samples included in Norfolk Island hybrid work as an example
$GATK_JAR --java-options "-Xmx4g" CombineGVCFs \
  -R reference.fasta \
  --variant Renunion_Zborbonics_WGS.15-179.g.vcf.gz \
  --variant N101.g.vcf.gz \
  --variant N102.g.vcf.gz \
  --variant N103.g.vcf.gz \
  --variant N104.g.vcf.gz \
  --variant N108.g.vcf.gz \
  --variant N110.g.vcf.gz \
  --variant N111.g.vcf.gz \
  --variant N114.g.vcf.gz \
  --variant N121.g.vcf.gz \
  --variant N122.g.vcf.gz \
  --variant N123.g.vcf.gz \
  --variant N126.g.vcf.gz \
  --variant N128.g.vcf.gz \
  --variant N129.g.vcf.gz \
  --variant N130.g.vcf.gz \
  --variant N208.g.vcf.gz \
  --variant N209.g.vcf.gz \
  --variant N212.g.vcf.gz \
  --variant NOR124.g.vcf.gz \
  --variant NOR125.g.vcf.gz \
  --variant NOR132.g.vcf.gz \
  --variant NOR134.g.vcf.gz \
  --variant NOR140.g.vcf.gz \
  --variant NOR141.g.vcf.gz \
  --variant NOR142.g.vcf.gz \
  --variant PN10.g.vcf.gz \
  --variant PN11.g.vcf.gz \
  --variant PN12.g.vcf.gz \
  --variant PN13.g.vcf.gz \
  --variant PN15.g.vcf.gz \
  --variant PN16.g.vcf.gz \
  --variant PN17.g.vcf.gz \
  --variant PN1.g.vcf.gz \
  --variant PN2.g.vcf.gz \
  --variant PN3.g.vcf.gz \
  --variant PN6.g.vcf.gz \
  --variant PN7.g.vcf.gz \
  -O NorfolkIsland_Hyrbidization.g.vcf.gz
