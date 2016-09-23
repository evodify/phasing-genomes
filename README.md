# Genome phasing
This is a pipeline to phase a genome into two homeologous subgenomes. This phasing approach requires parental genomes.

The phasing is performed in two stages:

1. Generate phased haplotype blocks using HapCUT (https://github.com/vibansal/hapcut)  
2. Merge haplotype block into two continuous strings.

Preparing reference parental genomes performed as an intermediate step. 

The pipeline is provided in the file HapCUT_phasing.sh

Required files:  

1) VCF file with samples to phase  
2) BAM files for samples to phase  
3) VCF file with parental genomes  
4) FASTA file with a reference genome  
