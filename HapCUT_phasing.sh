#!/bin/sh

##############################
# Select a sample from a VCF #
##############################

java -Xmx8g -jar GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R reference.fa \
  -V multiple_sample.vcf \
  -sn sample1 \
  -o sample1.vcf

##############
# Run HapCUT #
##############

# phase SNPs only

extractHAIRS --VCF sample1.vcf --bam sample1.bam --maxmem 128000 --mbq 20 --mmq 30 --PEonly 1 > sample1.fragment_matrix
HAPCUT --fragments sample1.fragment_matrix --VCF sample1.vcf --output sample1.haplotype --maxiter 100 --maxmem 128000  > sample1.haplotype.log
# --maxiter 100  and --maxiter 1000 produced the same results in my test.

# In my data, phasing SNPs-only data produced haplotypes that were more consistent with parental reference genomes that phasing both SNPs and indels. 
# to include indels, add options [--ref reference.fa --indels 1] to extractHAIRS command above
