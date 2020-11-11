# Novogene_WGS

The following repisitory contains scripts used for processing of Novogene WGS data produced for Zosterops lateralis (and other avian species)

1. FastQC.sh - conduct intial quality check using FastQC
2. Filter_Reads.sh - filter raw sequencing reads to remove low quality base calls and remove adapter content
3. FastQC2.sh - conduct secondary quality check using FastQC
4. Map_sample_reads.sh - map filtered reads to reference genome (outputs multiple bam files per individual)
5. Create_Single_Sample_Bam.sh - merge multiple bam files to generate a single bam for each individual
6. Sample_Genotyping.sh - conduct genotyping using GATK haplotype caller
7. Merge_gVCFs.sh - merge individual sample gVCFs into a multi-sample intermediate gVCF file

