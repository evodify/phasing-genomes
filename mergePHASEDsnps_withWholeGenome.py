#!/usr/bin/env python2
#
# The script merged phased SNPs dataset with whole genome homozygous sites with missing data correction.
# 
# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
#
# python mergePHASEDsnps_withWholeGenome.py -p {phased-data} -g {GT-data} -o {output} -Np 0.15
#
#*******************************************************************************************************

import argparse, numpy, sys

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

def hetToNs(g):
  if all([i in ['A','T','G','C','N'] for i in g]):
    gtN = numpy.random.binomial(1, prN, 1)
    gtl = [[[i, i], ['N', 'N']][int(gtN)] for i in g]
    gt = [i for e in gtl for i in e]
  else:
    gt = ['N', 'N']*len(g)
  return gt

phasedFile = open(args.phased, 'r')
phasedHeader = phasedFile.readline()
headerS = phasedHeader.split()
fileoutput = open(args.output, 'w')
fileoutput.write("%s\tREF\t%s\n" % ('\t'.join(str(e) for e in headerS[0:2]), '\t'.join(str(e) for e in headerS[2:])))

print('Merging ...')

phasedL = phasedFile.readline().split()
stopChr  = phasedL[0].split('_')[1]
stopPos = phasedL[1]
phasedGT = phasedL[2:]
GTfile = open(args.genome, 'r')
GTheader = GTfile.readline()

counter = 0
for line in GTfile:
  words = line.split()
  GTchr = words[0].split('_')[1]
  GTpos = words[1]
  chr_pos_ref = words[0:3]
  GTgt = words[3:]
  if stopChr == 'END':
    GTgtN = hetToNs(GTgt)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos):
    GTgtN = hetToNs(GTgt)
  elif (GTchr == stopChr and GTpos == stopPos):
    GTgtN = phasedGT
    phasedL = phasedFile.readline().split()
    if phasedL != []:
      stopChr  = phasedL[0].split('_')[1]
      stopPos = phasedL[1]
      phasedGT = phasedL[2:]
    else:
      stopChr = 'END'
      stopPos = []
      phasedGT = []
  else:
    print('ERROR')
    print('GT:', line)
    print('REF:', stopChr, stopPos)
  GTp = '\t'.join(str(e) for e in GTgtN)
  chr_pos_refP = '\t'.join(str(e) for e in chr_pos_ref)
  fileoutput.write("%s\t%s\n" % (chr_pos_refP, GTp))
  
  counter += 1
  if counter % 1000000 == 0:
    print str(counter), "lines processed"
  else:
    continue
    
GTfile.close()
phasedFile.close()
fileoutput.close()
print('Done!')
