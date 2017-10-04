#!/usr/bin/env python2
"""
This script rearrange haplotype blocks produced by HapCUT into two continuous strings. 

HapCUT_blocks.file:
BLOCK: offset: 3 len: 11 phased: 9 SPAN: 114 MECscore 19.92 fragments 43
3 0 1 scaffold_1  619 T TA  0/1:16,11:27:99:401,0,540 10,12:-53.7,-8.1,-27.1:-19.0:0.0
5 1 0 scaffold_1  655 TTA T 0/1:21,15:36:99:554,0,808 0,19:-78.9,-19.6,-4.1:15.5:4.0:FV
6 1 0 scaffold_1  656 TA  T 0/1:15,21:36:99:0|1:656_TA_T:805,0,586  7,0:-32.3,-32.4,-46.2:0.2:0.0
7 0 1 scaffold_1  658 TA  T 0/1:16,21:37:99:0|1:656_TA_T:805,0,586  15,0:-18.5,-32.5,-48.5:14.0:6.9:FV
8 0 1 scaffold_1  679 C A 0/1:15,21:36:99:793,0,519 15,23:-141.4,-60.0,-109.8:-49.8:1.0
9 0 1 scaffold_1  704 T C 0/1:19,8:27:99:274,0,739  21,8:-61.6,-63.7,-110.9:2.0:9.0
11  0 1 scaffold_1  722 C T 0/1:9,11:20:99:405,0,322  10,12:-91.3,-55.6,-84.2:-28.5:2.0
12  0 1 scaffold_1  727 CCATTAGAT C 0/1:8,8:16:99:310,0,312 5,0:-0.0,-0.0,-10.0:0.0:0.0
13  1 0 scaffold_1  733 G T,* 0/1:8,8,0:16:99:294,0,310,318,334,652 0,9:-34.4,-0.1,-0.1:0.0:0.0
******** 
BLOCK: offset: 16 len: 17 phased: 10 SPAN: 282 MECscore 5.98 fragments 10
16  0 1 scaffold_1  1194  A T 0/1:3,3:6:99:0|1:1194_A_T:117,0,117 3,3:-15.3,-7.8,-15.3:-7.5:0.0
17  0 1 scaffold_1  1204  G A 0/1:4,3:7:99:0|1:1194_A_T:114,0,159 4,3:-17.9,-13.0,-20.5:-4.9:1.0
18  1 0 scaffold_1  1210  A G 0/1:5,2:7:37:37,0,203 5,2:-15.6,-13.0,-23.4:-2.6:1.0
19  0 1 scaffold_1  1214  C T 0/1:4,4:8:99:0|1:1194_A_T:142,0,156 4,4:-23.1,-15.6,-23.1:-7.5:1.0
20  1 0 scaffold_1  1234  A T 0/1:6,2:8:37:37,0,161 6,2:-18.2,-15.6,-28.6:-2.6:1.0
25  0 1 scaffold_1  1297  C T 0/1:2,2:4:49:0|1:1297_C_T:78,0,49 3,2:-11.7,-10.4,-13.0:-1.3:2.0
26  0 1 scaffold_1  1301  C G 0/1:0,2:2:21:0|1:1297_C_T:81,0,21 0,3:-17.0,-7.8,-7.8:0.0:0.0
31  1 0 scaffold_1  1475  CAAAGAT C,* 0/1:9,1,0:10:1:0|1:1473_AT_A:1,0,332,28,335,363 2,0:-2.6,-2.6,-6.6:0.0:0.0
32  1 0 scaffold_1  1476  AAAGATG A,* 0/1:4,6,0:10:99:0|1:1474_TC_T:154,0,107,167,125,292 2,0:-2.6,-2.6,-6.6:0.0:0.0

parental_genomes_reference.file
# CHROM  POS parentA  parentB
scaffold_1  519 G GTTCGTCCCGTGTGAGTTCTCATTCAGCCATATTGTCAACTAACTCTGCGATATTCTTCACTC
scaffold_1  528 * A,G
scaffold_1  531 * C,T
scaffold_1  574 * C,T
scaffold_1  580 * A,T
scaffold_1  619 TA  T
scaffold_1  647 C A
scaffold_1  656 T TA
scaffold_1  658 T TA
scaffold_1  679 A C
scaffold_1  722 T C
scaffold_1  727 C CCATTAGAT
scaffold_1  1210  G A
scaffold_1  1234  T A
scaffold_1  1646  T A
scaffold_1  1651  A T
scaffold_1  1671  C A

output.file:
# CHROM  POS test.haplotype_A  test.haplotype_B
scaffold_1  619 TA  T
scaffold_1  655 N N
scaffold_1  656 TA  T
scaffold_1  658 N N
scaffold_1  679 A C
scaffold_1  704 C T
scaffold_1  722 T C
scaffold_1  727 C CCATTAGAT
scaffold_1  733 G T
scaffold_1  1194  A T
scaffold_1  1204  G A
scaffold_1  1210  G A
scaffold_1  1214  C T
scaffold_1  1234  T A
scaffold_1  1297  C T
scaffold_1  1301  C G
scaffold_1  1475  C CAAAGAT
scaffold_1  1476  A AAAGATG


command:
python assign_HapCUT_blocks.py -i HapCUT_blocks.file -r parental_genomes_reference.file -o output.file

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

# for options
import sys
import argparse
# for graphics
import matplotlib.pyplot as plt
import numpy as np

############################ options ##############################

class CommandLineParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = CommandLineParser()
parser.add_argument('-i', '--input_to_phase', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-r', '--phasing_reference', help = 'file containing list of alleles for parental species', type=str, required=True)
args = parser.parse_args()

############################ functions ###########################

def lists_overlap(gt, lst):
  ref = lst.split(",")
  return any(i == gt for i in ref)

def all_same(items):
    return all(x == items[0] for x in items)

def flat_list(lst):
  return [item for sublist in [i.split(",") for i in lst] for item in sublist]

def phase_state(GT, RefGT):
  overlapSame = [bool(lists_overlap(GT[0], RefGT[0])), bool(lists_overlap(GT[1], RefGT[1]))]
  overlapReverse = [bool(lists_overlap(GT[0], RefGT[1])), bool(lists_overlap(GT[1], RefGT[0]))]
  if overlapSame[0] and overlapSame[1]:
    return 'same'
  elif overlapReverse[0] and overlapReverse[1]:
    return 'reverse'

def phase_blocks(posBlock, GTblock, RefBlock, FlagB):
  """ rearrange HapCut blocks according to the order in parental reference file  """
  blockSameCount = 0
  blockReverseCount = 0
  GTblockPhase = []
  GTblockReturn = []

  for i in range(len(GTblock)):
    GT = GTblock[i]
    RefGT = RefBlock[i]
    if FlagB[i] == "FV": # uncertain variants are set to N
      GTblock[i] = ['N', 'N']
    else:  # find and count cases when phased genotype is consistent/inconsistent with parental genotypes
      GTphase = phase_state(GT, RefGT)
      if GTphase == 'same':
        blockSameCount += 1
        GTblockPhase.append('same')
      elif GTphase == 'reverse':
        blockReverseCount += 1
        GTblockPhase.append('reverse')

  # find prevalent phase 
  if all_same(GTblockPhase) and (len(GTblockPhase) >= 2):  # absolutely consistent with parental genotypes
    if GTblockPhase[0] == ['same']:
      RSratio = 1.0
    else:
      RSratio = 0.0
    RSratio = 0.0
  elif GTblockPhase == []:  # phase unknown
    RSratio = 'NA'
  else:
    RSratio = float(blockSameCount)/float(blockSameCount+blockReverseCount) # proportion of 'same' phasing state in block strings

  # define the block phase and produce output
  if (RSratio == 'NA') or (RSratio < 0.90 and RSratio > 0.10): # discard block that have > 90% of inconsistency with parental reference genotypes, or 
    for j in range(len(GTblock)):
      posBlockPrint = posBlock[j]
      GTblockPrint1 = 'N'
      GTblockPrint2 = 'N'
      GTblockReturn.append([posBlockPrint[0], posBlockPrint[1], GTblockPrint1, GTblockPrint2])
  else: # phase according to the prevalent state
    # find prevalent state
    phaseStateNumber = max(map(GTblockPhase.count, GTblockPhase))
    GTblockDefinedPahse = list(set( i for i in GTblockPhase if GTblockPhase.count(i) == phaseStateNumber ))
    if len(GTblockDefinedPahse) == 1:  # check if one state is prevalent
      if GTblockDefinedPahse == ['same']:
        phaseState = [0,1]
      else:
        phaseState = [1,0]
      for j in range(len(GTblock)):
        GT = GTblock[j]
        posBlockPrint = posBlock[j]
        GTblockPrint1 = GT[phaseState[0]]
        GTblockPrint2 = GT[phaseState[1]]
        GTblockReturn.append([posBlockPrint[0], posBlockPrint[1], GTblockPrint1, GTblockPrint2])
    else: # if there is conflict in phasing state, set to Ns. It usually applies for blocks with less then 10 position overlaps with parental reference.
      for j in range(len(GTblock)):
        posBlockPrint = posBlock[j]
        GTblockPrint1 = 'N'
        GTblockPrint2 = 'N'
        GTblockReturn.append([posBlockPrint[0], posBlockPrint[1], GTblockPrint1, GTblockPrint2])
        phaseState = [0,1]

  return(GTblockReturn, RSratio)

def write_phased(posB, GTb, RefB, Flag, R):
  """ preform phasing and write output"""
  output.write("********\n")  # uncomment this line to separate blocks
  phasedBlockRatio = phase_blocks(posB, GTb, RefB, Flag)
  phasedBlock = phasedBlockRatio[0]
  if not (phasedBlockRatio[1] == 0 or phasedBlockRatio[1] == 1 or phasedBlockRatio[1]=='NA'): # ignore ratio of 0, 1, NA
    R.append(phasedBlockRatio[1])
  for block in phasedBlock:
    BlockPrint = '\t'.join(str(e) for e in block)
    output.write("%s\n" % BlockPrint)

############################ script ##############################

RefFileName = args.phasing_reference
RefFile = open(RefFileName, "r")
refHead = RefFile.readline()
refWords = RefFile.readline().split()
refChr = int(refWords[0].split("_")[1])
refPos = int(refWords[1])

counter = 0

print('Creating the output file...')
sampleID = args.input_to_phase
output = open(args.output, 'w')
output.write("#CHROM\tPOS\t%s_A\t%s_B\n" % (sampleID, sampleID)) # make a header

Ratio = [] # ratio between phasing states ['reverse', 'same']

with open(args.input_to_phase) as datafile:
  for line in datafile:
    words = line.split()

    # read blocks, phase and write output
    if words[0] == "********": # if the end of a block is reached, phase the block and write output
      write_phased(posBlock, GTblock, RefBlock, FlagBlock, Ratio)

    # read a block and corresponding reference info into memory
    elif words[0] == "BLOCK:":  # reset all lists at the beginning of each block
      GTblock = []
      posBlock = []
      RefBlock = []
      FlagBlock = []
    else:  # read a block
      phase = words[1:3]
      chrPos = words[3:5]
      genotypes = flat_list(words[5:7])
      GT = [genotypes[int(phase[0])], genotypes[int(phase[1])]]
      # check whether a variant marked as not real or not heterozygous (FV flag by HapCUT)
      mec = words[8].split(":")
      if mec[-1] == 'FV': 
        flag = 'FV'
      else:
        flag = 'no' 

      # find overlap with parental reference genotypes
      gtChr = int(chrPos[0].split("_")[1])
      gtPos = int(chrPos[1])
      # read reference file until the overlap with input data is found
      while refChr < gtChr or refChr == gtChr and refPos < gtPos:
        refWords = RefFile.readline().split()
        if not refWords:
          break
        else:
          refChr = int(refWords[0].split("_")[1])
          refPos = int(refWords[1])
      if gtChr == refChr and gtPos == refPos:  # set reference genotypes 
        RefGT = refWords[2:]
      else:  # if not found, set reference to Ns
        RefGT = ['N', 'N']

      # append the genotypes, positions, reference and hapCUT flags
      GTblock.append(GT)
      posBlock.append(chrPos)
      RefBlock.append(RefGT)
      FlagBlock.append(flag)

    # track the progress:
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

  # Phase the last block
  write_phased(posBlock, GTblock, RefBlock, FlagBlock, Ratio)

  output.close()
datafile.close()
RefFile.close()

# plot the Ratio distribution
plt.hist(Ratio, color="grey", bins = 100)
plt.xticks(np.arange(0,1, 0.1))
plt.xlabel('Propostion of phasing states: same/(same+reverse)')
plt.ylabel("Number of blocks")
plt.savefig(args.output+'.png', dpi=90)

print "Done!"
