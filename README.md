# Genome phasing
This is a pipeline to phase a genome into two homeologous subgenomes. Parental genomes are required.

The phasing is performed in two stages:

1. Generate phased haplotype blocks using [HapCUT](https://github.com/vibansal/hapcut)
2. Merge haplotype blocks into two continuous strings.

Preparing reference parental genomes performed as an intermediate step. 

Required files:  

1) VCF file with samples to phase  
2) BAM files for samples to phase  
3) VCF file with parental genomes  
4) FASTA file with a reference genome  

## Generate phased haplotype blocks using HapCUT

### Select a sample from a multisample VCF
Can be performed with [GATK](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php) :
```
java -Xmx8g -jar GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R reference.fa \
  -V multiple_sample.vcf \
  -sn sample1 \
  -o sample1.vcf
```

### Run HapCUT

In my experience, phasing SNPs-only data produced haplotypes that were more consistent with parental reference genomes than phasing both SNPs and indels. If indels information is not required, I recommend to use only SNPs data.

```
extractHAIRS --VCF sample1.vcf --bam sample1.bam --maxmem 128000 --mbq 20 --mmq 30 --PEonly 1 > sample1.fragment_matrix
HAPCUT --fragments sample1.fragment_matrix --VCF sample1.vcf --output sample1.haplotype --maxiter 100 --maxmem 128000  > sample1.haplotype.log
```
To include indels, add to `extractHAIRS` command above the following options: `--ref reference.fa --indels 1`

## Generate parental reference genomes

### Convert VCF to tab-deilimited table

Performed with [GATK](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php) :
```
java -Xmx8g -jar GenomeAnalysisTK.jar \
 -T VariantsToTable \
 -R reference.fa \
 -V reference_genomes_GT.vcf \
 -F CHROM -F POS -F REF -F ALT -GF GT \
 -o reference_genomes_GT.table
```
 
### Make reference file
```
python createREFgenomesForPhasing.py -i reference_genomes_GT.table -o reference_genomes_REF.tab -s1 parent1_1,parent1_2 -s2 parent2_1,parent2_2  -m 0.25
```
### Keep fixed differences only
```
python filterREFgenomeFixedOnly.py -i reference_genomes_REF.tab -o reference_genomes_REF.fixedOnly.tab
```
## Merge haplotype blocks
```
python merge_HapCUT_blocks.py -i sample1.haplotype -r reference_genomes_REF.fixedOnly.tab -o sample1.haplotype.PHASED
grep -v "\*\*\*" sample1.haplotype.PHASED > sample1.haplotype.PHASED.tab
```
