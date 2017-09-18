#!/usr/bin/env python2
'''
The script merged phased SNPs dataset with whole genome homozygous sites. Because some heterozygouse sites are not phased and set to Ns, missing data shoudl be also introduced to homozygous sites. This missing data correction is performed with -Np option.

phased-file:
CHROM   POS Sample1_A   Sample1_B   Sample2_A   Sample2_B
scaffold_1  2   G   T   G   T
scaffold_1  5   A   T   A   T
scaffold_1  7   G   A   G   A
scaffold_1  8   G   A   G   A
scaffold_1  11  C   T   C   T

GT-file can be one- or two-character coded.
two-character coded GT-file:
CHROM   POS Sample1 Sample2
scaffold_1  1   A/A A/A
scaffold_1  2   A/T T/T
scaffold_1  3   T/T T/T
scaffold_1  4   A/A A/A
scaffold_1  5   G/G G/G
scaffold_1  6   C/C C/C
scaffold_1  7   G/A G/A
scaffold_1  8   G/A G/A
scaffold_1  9   G/G G/G
scaffold_1  10  T/T T/T
scaffold_1  11  C/T C/T
scaffold_1  12  A/A A/T
scaffold_1  13  G/G G/G
scaffold_1  14  T/T T/T
scaffold_1  15  C/C C/C

one-character coded GT-file:
CHROM   POS     Sample1 Sample2
scaffold_1      1       A       A
scaffold_1      2       W       T
scaffold_1      3       T       T
scaffold_1      4       A       A
scaffold_1      5       G       G
scaffold_1      6       C       C
scaffold_1      7       R       R
scaffold_1      8       R       R
scaffold_1      9       G       G
scaffold_1      10      T       T
scaffold_1      11      Y       Y
scaffold_1      12      A       W
scaffold_1      13      G       G
scaffold_1      14      T       T
scaffold_1      15      C       C

output-file:
CHROM   POS Sample1_A   Sample1_B   Sample2_A   Sample2_B
scaffold_1  1   A   A   A   A
scaffold_1  2   G   T   G   T
scaffold_1  3   T   T   T   T
scaffold_1  4   A   A   A   A
scaffold_1  5   A   T   A   T
scaffold_1  6   C   C   C   C
scaffold_1  7   G   A   G   A
scaffold_1  8   G   A   G   A
scaffold_1  9   G   G   G   G
scaffold_1  10  T   T   T   T
scaffold_1  11  C   T   C   T
scaffold_1  12  N   N   N   N
scaffold_1  13  G   G   G   G
scaffold_1  14  T   T   T   T
scaffold_1  15  C   C   C   C


contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

python mergePHASEDsnps_withWholeGenome.py -p phased-file -g GT-file -o output-file -Np 0.15

'''

############################# modules #############################

import sys, argparse  # for options
import numpy  # to assign Ns with probability

############################ options ##############################
class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)
      
parser = MyParser()
parser.add_argument('-p', '--phased', help = 'name of the phased data file', type=str, required=True)
parser.add_argument('-g', '--genome', help = 'name of the GT table file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-Np', '--probability_of_N', help = 'probability of N in homozygotes', type=float, required=True)
args = parser.parse_args()

prN = args.probability_of_N

############################ functions ##############################

def hetToNs(g, pN):
  ''' all unphased heterozygotes are set to Ns, homozygotes are set to Ns with probability pN '''
  if all([i in ['A','T','G','C','N'] for i in g]):
    gtN = numpy.random.binomial(1, pN, 1)
    gtl = [[[i, i], ['N', 'N']][int(gtN)] for i in g]
    gt = [i for e in gtl for i in e]
  else:
    # sites with heterozygotes from GT-data are not printed
    # if printed by mistake, such sites are set to all-Ns
    gt = ['N', 'N']*len(g)
  return gt

def is_polymorphic(sampWords):
  ''' check if the set of genotypes is polymorphic '''
  # fist skip missing data
  noNsGT = []
  for i in (sampWords):
    if i != 'N':
      noNsGT.append(i)
  # check if there is polymorphism:
  return any(x in 'RYMKSW' or x != noNsGT[0] for x in noNsGT)

def process_unphased(g, pN):
  """ replace heterozygouts with Ns and skip polymorphic sites in unphased data"""
  gN = hetToNs(g, pN)
  if is_polymorphic(gN):
    gN = ['N', 'N']*len(g)
  return gN

############################ script ##############################

phasedFile = open(args.phased, 'r')
phasedHeader = phasedFile.readline()
headerS = phasedHeader.split()
fileoutput = open(args.output, 'w')
fileoutput.write("%s\t%s\n" % ('\t'.join(str(e) for e in headerS[0:2]), '\t'.join(str(e) for e in headerS[2:])))

print('Merging ...')

phasedL = phasedFile.readline().split()
stopChr  = phasedL[0].split('_')[1]
stopPos = phasedL[1]
phasedGT = phasedL[2:]
GTfile = open(args.genome, 'r')
GTheader = GTfile.readline().split()

# check if number of samples in two files overlap
if len(phasedGT) != len(GTheader[2:])*2:
  raise Exception('Number of samples differs in two files')

counter = 0

for line in GTfile:
  words = line.split()
  GTchr = words[0].split('_')[1]
  GTpos = words[1]
  chr_pos = words[0:2]

  # define is GT.table is one- or two-character code
  if '/' in words[3]:
    GTgt = []
    # set to one-character code
    GT = [i.split("/") for i in words[2:]]
    for i in GT:
      if i[0] != i[1]:
        GTgt.append('H') # all heterozygotes are set to H
      else:
        GTgt.append(i[0])
  else:
    GTgt = words[2:]

  # merge files
  if stopChr == 'END':
    # if the end of phased data reached, introduce Ns to the last block 
    GTgtN = process_unphased(GTgt, prN)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos):
    # read line of GT.table while it is on the position before the phased position
    GTgtN = process_unphased(GTgt, prN)
  elif int(GTchr) < int(stopChr):
    # read line of GT.table while it is on the chromosome before the phased chromosome
    GTgtN = process_unphased(GTgt, prN)
  elif (GTchr == stopChr and GTpos == stopPos):
    # if at the same position, read next line of phased data
    GTgtN = phasedGT
    phasedL = phasedFile.readline().split()
    if phasedL != []: # if phased data present, read phased line
      stopChr  = phasedL[0].split('_')[1]
      stopPos = phasedL[1]
      phasedGT = phasedL[2:]
    else: # if the end of phased data reached
      stopChr = 'END'
      stopPos = []
      phasedGT = []
  else:
    print line
    # print Error with positions where if happened if there is an unpredicted condition
    print('ERROR')
    print('GT:', GTchr, GTpos)
    print('REF:', stopChr, stopPos)
  
  # write output
  GTp = '\t'.join(str(e) for e in GTgtN)
  chr_pos_refP = '\t'.join(str(e) for e in chr_pos)
  fileoutput.write("%s\t%s\n" % (chr_pos_refP, GTp))

  # track the progress:
  counter += 1
  if counter % 1000000 == 0:
    print str(counter), "lines processed"

GTfile.close()
phasedFile.close()
fileoutput.close()

print('Done!')
