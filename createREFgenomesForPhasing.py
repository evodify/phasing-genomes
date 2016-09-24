#!/usr/bin/python2
#
# This script makes a references of two parental genomes by summarizing alleles.
#
# datafile:
# #CHROM POS REF ALT ind1 ind2 ind3 ind4
# scaffold_1  1 A T ./. ./. A/A A/A
# scaffold_1  12 G A,* A/A A/A G/G G/G
# scaffold_1  54 C T,* */* */* C/T C/C
# scaffold_1  81 T TA  TA/TA TA/TA T/T T/T
#
# outputfile:
# #CHROM  POS Species1  Species2
# scaffold_1  12 A G
# scaffold_1  54 * C,T
# scaffold_1  81 TA T
#
# command:
# python createREFgenomesForPhasing.py -i <datafile> -o <outputfile> -s1 <ind1,ind2> -s2 <ind3,ind4> -m 0.1
#
# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
#
#*******************************************************************************************************
import argparse, re, collections, sys

def removeNgap(items): # to remove N and gaps
  items2 = [x for x in items if not (x == './.')]
  return items2
  
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

miss = args.missingness
print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  ind_1 = [str(j) for j in re.findall(r'[^,\s]+', args.samples1)]
  ind_sp1 = []
  for i in ind_1:
    indnumber = header_words.index(i)
    ind_sp1.append(indnumber)
  
  ind_2 = [str(j) for j in re.findall(r'[^,\s]+', args.samples2)]
  ind_sp2 = []
  for i in ind_2:
    indnumber = header_words.index(i)
    ind_sp2.append(indnumber)

  print('Creating the output file...')
  sp12_output = open(args.output, 'w')
  sp12_output.write("#CHROM\tPOS\tSpecies1\tSpecies2\n")

  print('Creating reference genomes...')
  counter = 0
  conterProcessed = 0
  conterMissing = 0

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    rowName = '\t'.join(str(e) for e in chr_pos)
    sp1_alleles = []
    for el in ind_sp1:
      sp1_alleles.append(words[el])

    sp2_alleles = []
    for el in ind_sp2:
      sp2_alleles.append(words[el])

    sp1_alleles_noN = removeNgap(sp1_alleles)
    sp2_alleles_noN = removeNgap(sp2_alleles)
    sp1_alleles_noNspl = [item for sublist in [i.split("/") for i in sp1_alleles_noN] for item in sublist]
    sp2_alleles_noNspl = [item for sublist in [i.split("/") for i in sp2_alleles_noN] for item in sublist]    
    if float(len(sp1_alleles_noN)) >= float(len(sp1_alleles))*miss and float(len(sp2_alleles_noN)) >= float(len(sp2_alleles))*miss:
      conterProcessed += 1
      sp1_alleles_noNsplset = set(sp1_alleles_noNspl)
      sp2_alleles_noNsplset = set(sp2_alleles_noNspl)
      sp1AltPrint = ','.join(str(al) for al in sp1_alleles_noNsplset)
      sp2AltPrint = ','.join(str(al) for al in sp2_alleles_noNsplset)
      sp12_output.write("%s\t%s\t%s\n" % (rowName, sp1AltPrint, sp2AltPrint))
    else:
      conterMissing += 1
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"
    else:
      continue
  print('Done!\n')
  print str(conterProcessed), "sites processed"
  print str(conterMissing), "sites with too many missing genotypes in both species"
datafile.close()
sp12_output.close()

