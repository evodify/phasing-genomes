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

All python script contain description of input and output data format in the header of each file.
To see possible option, run python script with --help option:
`python script.py --help`

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

In my experience, phasing SNPs-only data produced haplotypes that were more consistent with parental reference genomes than phasing both SNPs and indels. If indels information is not required, I recommend using only SNPs data.

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
`multiple_sample.vcf` should also be converted to `multiple_sample_GT.table` using this approach.

### Make a reference file
```
python createREFgenomesForPhasing.py -i reference_genomes_GT.table -o reference_genomes_REF.tab -s1 parent1_1,parent1_2 -s2 parent2_1,parent2_2  -m 0.25
```
### Keep non-shared polymorphism only
```
python filterREFgenomeNonSharedOnly.py -i reference_genomes_REF.tab -o reference_genomes_REF.nonShared.tab
```
## Phase haplotype blocks
```
python assign_HapCUT_blocks.py -i sample1.haplotype -r reference_genomes_REF.nonShared.tab -o sample1.haplotype.PHASED
```
Chimeric blocks (the phasing state was supported by less than 90% of sites) were set to missing data.
For example, below is the distribution of phasing state. Blocks between 0.10 and 0.90 are considered chimeric in the script. If your distrubution is different (graphics is produced by the script), edit this code `RSratio < 0.90 and RSratio > 0.10` in the script.

![alt tag] (https://github.com/evodify/genome-phasing/blob/master/chimericBlocks.png)

### Merge heterozygous and homozygous sites

`sample1.haplotype.PHASED` contains only heterozygous sites of sample1. However, originally sample1 also contained homozygous sites that were polymorphic in other samples. These homozygous sites need to be returned.
Phasing introduces some amount of missing data. To keep balance between homozygous and heterozygous sites, the amount of introduced Ns need to be assessed and the same amount of Ns should be introduced to homozygous sites.

#### Estimate the missing data correction value.
```
# heterozygotsPhased:
for i in *.haplotype.PHASED; do grep -vwc "********" $i; done
# nonMissingPhased:
for i in *.haplotype.PHASED; do grep -vw "********" $i | awk '$3!=$4 {count++} END {print count}'; done
```

`introducedNs = 1 - nonMissingPhased / heterozygotsPhased`

`introducedNs` is used to define the missing data correction value (-Np). In my case, it was ~0.20.

#### Merge with introduction of Ns
```
python mergePhasedHeteroHomo_randomNs.py -p sample1.haplotype.PHASED -s sample-name -g multiple_sample_GT.table -o sample1.haplotype.PHASED.tab -Np 0.20
```

##### Count heterozygous and homozygous sites in original files
`n` should be replaced with number of samples + 3:
```
# homozygotsOriginal:
for i in {4..n}; do cut -f $i multiple_sample_GT.table | sed 's/\// /g;s/\./N/g' | awk '$1==$2 && $1!="N" && $2!="N" {count++} END {print count}'; done

# heterozygotsOriginal:
for i in {4..n}; do cut -f $i multiple_sample_GT.table | sed 's/\// /g;s/\./N/g' | awk '$1!=$2 {count++} END {print count}'; done
```

##### Count heterozygous and homozygous sites in merged files
```
# homozygotsPhased:
for i in *.haplotype.PHASED.tab; do awk '$3==$4 && $3!="N" && $4!="N" {count++} END {print count}' $i; done

# heterozygotsPhased:
for i in *.haplotype.PHASED.tab; do awk '$3!=$4 {count++} END {print count}' $i; done
```

##### Compare the levels of heterozygosity between original and phased data. 
Expectation: `heterozygotsOriginal / homozygotsOriginal = heterozygotsPhased / homozygotsPhased`


The value of `introducedNs` is used for a very rough correction. For the most precise correction, `introducedNs` should be lowered by some amount because extra Ns are introduced to homozygots in all-Ns blocks.
I recommend to run `mergePhasedHeteroHomo_randomNs.py` with `introducedNs` values in the -Np option and then lower it little by little, until the ratio heterozygots / heterozygots in the phased and non-phased data are as similar as possible.

### Merge all phased files togather
```
for i in *.haplotype.PHASED.tab; do cut -f 3,4 $i.col34; done
rm sample1.haplotype.PHASED.tab.col34
paste 12.4.haplotype.PHASED.tab *.col34 > all.haplotype.PHASED.tab
```

### Merge phased SNPs with a whole genome (optional)

```
python mergePHASEDsnps_withWholeGenome.py -p all.haplotype.PHASED.tab -g whole_genome_multiple_sample_GT.tab -o all.haplotype.PHASED.wholeGenome.tab -Np 0.16
```
Only homozygous sites from a whole genome will be used for merging. Unphased heterozygous sites will be set to Ns. Number of samples in `all.haplotype.PHASED.tab` and `whole_genome_multiple_sample_GT.tab` should be the same.
Again, missing data is a problem here. Phasing introduced some amount of Ns, so this needs to be taken into account during merging with a whole genome. MissingCorrectionValue (0.16) also need to be used here. 

**Note!** Check the ratio between polymorphic and non-polymorphic sites before and after phasing. It should be the same. If it is not, modify MissingCorrectionValue until you get the same proportion. Artificially changing polymorphic/non-polymorphic ration can biase results in some subsequent analyses.
