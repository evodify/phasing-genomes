# Genome phasing
This is a pipeline to phase a genome into two homeologous subgenomes. This phasing approach requires parental genomes.

The phasing is performed in two stages:

1) Generate phased haplotype blocks using HapCUT (https://github.com/vibansal/hapcut)  
2) Merge haplotype block into two continuous strings.

Preparing reference parental genomes performed as an intermediate step. 

The pipeline is provided in the file HapCUT_phasing.sh


