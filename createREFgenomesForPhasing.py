#!/usr/bin/env python2
'''
This script makes a references of two parental genomes by summarizing alleles.

inputfile:
#CHROM POS REF ALT ind1 ind2 ind3 ind4
scaffold_1  1 A T ./. ./. A/A A/A
scaffold_1  12 G A,* A/A A/A G/G G/G
scaffold_1  54 C T,* */* */* C/T C/C
scaffold_1  81 T TA  TA/TA TA/TA T/T T/T

outputfile:
#CHROM  POS Species1  Species2
scaffold_1  12 A G
scaffold_1  54 * C,T
scaffold_1  81 TA T

command:
python createREFgenomesForPhasing.py -i inputfile -o outputfile -s1 ind1,ind2 -s2 ind3,ind4 -m 0.1

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

import argparse, re, collections, sys

############################ options ##############################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s1', '--samples1', help = 'column names of the species A (comma delimited)', type=str, required=True)
parser.add_argument('-s2', '--samples2', help = 'column names of the species B (comma delimited)', type=str, required=True)
parser.add_argument('-m', '--missingness', help = 'allowed fraction of missing data per site per species', type=float, required=True)
args = parser.parse_args()

############################ functions ###########################

def index_sample(header, samplList):
  ''' index position of samples in file '''
  ind_names = [str(j) for j in re.findall(r'[^,\s]+', samplList)]
  ind_index = []
  for i in ind_names:
    indnumber = header.index(i)
    ind_index.append(indnumber)
  return ind_index

def extract_genotypes(genotypes, ind_index):
  ''' extract genotypes using samples index '''
  alleles = []
  for g in ind_index:
    alleles.append(genotypes[g])
  # remove missing data:
  allelesNoNs = [x for x in alleles if not (x == './.')]
  # split alleles:
  allelesNoNsSplit = [i for sublist in [i.split("/") for i in allelesNoNs] for i in sublist]
  return allelesNoNsSplit

############################ script ##############################

miss = args.missingness
counter = 0
counterProcessed = 0
counterMissing = 0

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()
  
  # index samples of the two species
  ind_sp1 = index_sample(header_words, args.samples1)
  ind_sp2 = index_sample(header_words, args.samples2)

  print('Creating the output file...')
  sp12_output = open(args.output, 'w')
  sp12_output.write("#CHROM\tPOS\tSpecies1\tSpecies2\n")

  print('Creating reference genomes...')
  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    rowName = '\t'.join(str(e) for e in chr_pos) # make chr pos row names
    
    # extract genotypes using samples index
    sp1_alleles = extract_genotypes(words, ind_sp1)
    sp2_alleles = extract_genotypes(words, ind_sp2)

    if float(len(sp1_alleles)) >= float(len(sp1_alleles))*miss and float(len(sp2_alleles)) >= float(len(sp2_alleles))*miss:
      counterProcessed += 1
      sp1_allelesSet = set(sp1_alleles)
      sp2_allelesSet = set(sp2_alleles)
      sp1AltPrint = ','.join(str(al) for al in sp1_allelesSet)
      sp2AltPrint = ','.join(str(al) for al in sp2_allelesSet)
      sp12_output.write("%s\t%s\t%s\n" % (rowName, sp1AltPrint, sp2AltPrint))
    else:
      counterMissing += 1
    
    # track the progress:
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"
      
datafile.close()
sp12_output.close()    
print('Done!\n')

# print statistics on the screen
print str(counterProcessed), "sites processed"
print str(counterMissing), "sites with too many missing genotypes in both species\n"
