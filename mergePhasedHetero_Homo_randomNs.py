#!/usr/bin/env python2
'''
The script merged phased heterozygous sites with homozygous sites that were skipped during phasing by HapCUT.

phased-data:
CHROM  POS Sample1_A  Sample1_B
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

#GT.table can be one- or two-character coded.

# two-character coded GT.table:

CHROM	POS	Sample1	Sample2	Sample3	Sample4
scaffold_1	508	T/T	T/T	T/T	A/A
scaffold_1	593	T/T	T/T	T/A	T/T
scaffold_1	619	TA/T	TA/T	TA/T	TA/T
scaffold_1	635	A/A	A/A	A/T	A/A
scaffold_1	655	N/N	N/N	A/T	A/A
scaffold_1	656	TA/T	TA/T	TA/T	TA/T
scaffold_1	658	N/N	N/N	A/A	T/A
scaffold_1	679	A/C	A/C	A/C	A/C
scaffold_1	700	C/C	C/C	C/C	C/T
scaffold_1	704	C/T	C/T	C/T	C/T
scaffold_1	722	T/C	T/C	T/C	T/C
scaffold_1	727	C/CCATTAGAT	C/CCATTAGAT	C/CCATTAGAT	C/CCATTAGAT
scaffold_1	733	G/T	G/T	G/T	G/T
scaffold_1	900	T/G	G/G	G/G	G/G
scaffold_1	1000	A/T	T/T	T/T	T/T
scaffold_1	1194	A/T	A/T	A/T	A/T
scaffold_1	1204	G/A	G/A	G/A	G/A
scaffold_1	1210	G/A	G/A	G/A	G/A
scaffold_1	1214	C/T	C/T	C/T	C/T
scaffold_1	1234	T/A	T/A	T/A	T/A
scaffold_1	1297	C/T	C/T	C/T	C/T
scaffold_1	1301	C/G	C/G	C/G	C/G
scaffold_1	1475	C/CAAAGAT	C/CAAAGAT	C/CAAAGAT	C/CAAAGAT
scaffold_1	1476	A/AAAGATG	A/AAAGATG	A/AAAGATG	A/AAAGATG

# one-character coded GT.table (NOTE! Skips indels):

CHROM   POS Sample1 Sample2 Sample3 Sample4
scaffold_1  508 T   T   T   A
scaffold_1  593 T   T   W   T
scaffold_1  619 N   N   N   N
scaffold_1  635 A   A   W   A
scaffold_1  655 N   N   W   A
scaffold_1  656 N   N   N   N
scaffold_1  658 N   N   A   W
scaffold_1  679 M   M   M   M
scaffold_1  700 C   C   C   Y
scaffold_1  704 Y   Y   Y   Y
scaffold_1  722 Y   Y   Y   Y
scaffold_1  727 N   N   N   N
scaffold_1  733 K   K   K   K
scaffold_1  900 K   G   G   G
scaffold_1  1000    W   T   T   T
scaffold_1  1194    W   W   W   W
scaffold_1  1204    R   R   R   R
scaffold_1  1210    R   R   R   R
scaffold_1  1214    Y   Y   Y   Y
scaffold_1  1234    W   W   W   W
scaffold_1  1297    Y   Y   Y   Y
scaffold_1  1301    S   S   S   S
scaffold_1  1475    N   N   N   N
scaffold_1  1476    N   N   N   N

output:
CHROM  POS Sample1_A  Sample1_B
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

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

python mergePhasedHetero_Homo_randomNs.py -p phased-data -g GT.table -s sample-name -o output  -Np 0.15
'''

import argparse, sys, numpy

############################ options ##############################

class CommandLineParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = CommandLineParser()
parser.add_argument('-p', '--phased', help = 'name of the phased data file', type=str, required=True)
parser.add_argument('-g', '--gt_table', help = 'name of the GT table file', type=str, required=True)
parser.add_argument('-s', '--sample_name', help = 'name of the sample in the GT table file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-Np', '--probability_of_N', help = 'probability of N in homozygotes', type=float, required=True)
args = parser.parse_args()

prN = args.probability_of_N

############################ functions ###########################

def oneToTwoCharacter(g):
  ''' splits one-character code to two-characters coded genotypes '''
  if g in 'ATGC':
    gt = [g]*2
  elif g == 'R':
   gt = ['A','G']
  elif g == 'Y':
    gt = ['T','C']
  elif g == 'M':
    gt = ['A','C']
  elif g == 'K':
    gt = ['G','T']
  elif g == 'S':
    gt = ['G','C']
  elif g == 'W':
    gt = ['A','T']
  else:
    gt = ['N', 'N']
  return gt

def hetToNsBlock(g):
  ''' heterozygotes from GT.table that are within a Ns block are set to Ns '''
  if g[0] != g[1] or g == ['.', '.']:
    g = ['N', 'N']
  return g

def hetToNsBefore(g):
  ''' all heterozygotes and some homozygotes are set to Ns '''
  if g[0] != g[1] or g == ['.', '.']:
    g = ['N', 'N']
  else:
    gtN = numpy.random.binomial(1, prN, 1)
    g = [g, ['N', 'N']][int(gtN)]
  return g

def writeBlock(listToWrire, fileout):
  ''' writes an assembled block to a file '''
  if listToWrire:
    if all(x[1][0] == x[1][1] for x in listToWrire):
      for i in listToWrire:
        i[1] = ['N', 'N']
    #fileout.write("******* start block\n")  # for debug mode
    for el in range(len(listToWrire)):
      listToWrireP = '\t'.join(str(i) for sublist in listToWrire[el] for i in sublist)
      fileout.write("%s\n" % listToWrireP)
    #fileout.write("******* end block\n")  # for debug mode

def writeBefore(listToWrire, fileout):
  ''' writes into a file homozygotes that are between blocks '''
  if listToWrire:
    #fileout.write("******* start inter-block\n")  # for debug mode
    for el in range(len(listToWrire)):
      if listToWrire[el][1][0] != listToWrire[el][1][1] or listToWrire[el][1] == ['.', '.']:
        listToWrire[el][1] = ['N', 'N']
      listToWrireP = '\t'.join(str(i) for sublist in listToWrire[el] for i in sublist)
      fileout.write("%s\n" % listToWrireP)
    #fileout.write("******* end inter-block\n")  # for debug mode
    
def appendGT(coord, g, block):
  ''' appends genotypes from GT.table and sets heterozygotes to Ns on a fly'''
  if g[0] != g[1] or g == ['.', '.']:
    gt = ['N', 'N']
  else:
    gt = g
  block.append([coord, g])
  
def readPhasedLine(fileName, stat):
  ''' reads lines from the phased file with block separation '''
  line = fileName.readline().split()
  if line == []:
    pChr = 'NA'
    pPos = 'NA'
    pGT = 'NA'
    stat = 'before'
  elif line == ['********']:
    line = fileName.readline().split()
    pChr  = line[0].split('_')[1]
    pPos = line[1]
    pGT = line[2:4]
    if stat == 'before':
      stat = 'block'
    else:
      GTblock = []
      stat = 'before'
  else:
    pChr  = line[0].split('_')[1]
    pPos = line[1]
    pGT = line[2:4]
  return (pChr, pPos, pGT, stat)

def appendPhasedGT(GTblo, coord, pGT, pFile):
  ''' appends genotypes from a phased file'''
  stat = 'block'
  GTblo.append([coord, pGT])
  pL = readPhasedLine(pFile, stat)
  pChr  = pL[0]
  pPos = pL[1]
  pGT = pL[2]
  stat = pL[3]
  return (GTblo, pChr, pPos, pGT, stat)

def writeFragment(stat1, stat, GTbe, GTblo, fileout):
  ''' decides on whether a block or inter-block fragments and writes to the file '''
  if stat1 == 'before' and stat == 'block':
    writeBefore(GTbe, fileout)
    GTbe = []
    stat1 = stat
  elif stat == 'before' and stat1 == 'block':
    writeBlock(GTblo, fileout)
    GTblo = []
    stat1 = stat
  return (stat1, GTbe, GTblo)

############################ script ##############################

name = args.sample_name

phasedFile = open(args.phased, 'r')
phasedHeader = phasedFile.readline()
fileoutput = open(args.output, 'w')
fileoutput.write("%s" % phasedHeader)

print('Merging ...')

status = 'before'
status1 = 'before'
GTbefore = []
GTblock = []

phasedLine1 = phasedFile.readline().split()
phasedL = readPhasedLine(phasedFile, status)
stopChr  = phasedL[0]
stopPos = phasedL[1]
phasedGT = phasedL[2]

GTfile = open(args.gt_table, 'r')
GTheader = GTfile.readline().split()
indnumber = GTheader.index(name) # index sample position

counter = 0

for line in GTfile:
  words = line.split()
  GTchr = words[0].split('_')[1]
  GTpos = words[1]
  GTcoord = words[0:2]

  # define is GT.table is one- or two-character code
  if '/' in words[indnumber]:
    GTgt = words[indnumber].split('/')
  else:
    GTgt = oneToTwoCharacter(words[indnumber])

  if stopChr == 'NA':
    ''' when there is no phased data but GT.table still has some information. 
    Usually it is at the end of file '''
    GTblock = []
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos) and status == 'before':
    # append information from GT.table that is between phased blocks
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos) and status == 'block':
    # append information from GT.table that is within a phased block
    GTgtN = hetToNsBlock(GTgt)
    appendGT(GTcoord, GTgtN, GTblock)
  elif (GTchr == stopChr and GTpos == stopPos):
    # append information from a phased file
    eqLine = appendPhasedGT(GTblock, GTcoord, phasedGT, phasedFile)
    GTblock = eqLine[0]
    stopChr = eqLine[1]
    stopPos = eqLine[2]
    phasedGT = eqLine[3]
    status = eqLine[4]
    if GTblock:
      wF = writeFragment(status1, status, GTbefore, GTblock, fileoutput)
    status1 = wF[0]
    GTbefore = wF[1]
    GTblock = wF[2]
  elif int(GTchr) != int(stopChr):
    # write output if the end of a chromosome is reached
    status = 'before'
    writeBlock(GTblock, fileoutput)
    GTblock = []
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) > int(stopPos):  
    # read more lines from GT.table if its position exceed the phased data
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
    # read phased file untill the same position is reached
    while (int(GTchr) > int(stopChr)) or (int(GTchr) == int(stopChr) and int(GTpos) > int(stopPos)):
      phasedL = readPhasedLine(phasedFile, status)
      stopChr  = phasedL[0]
      stopPos = phasedL[1]
      phasedGT = phasedL[2]
    # write inter-block fragment
    wF = writeFragment(status1, status, GTbefore, GTblock, fileoutput)
    status1 = wF[0]
    GTbefore = wF[1]
    GTblock = wF[2]
  else: # print Error if there is an unpredicted condition
    raise Exception('ERROR: unexpected case, check input files or the script for errors')

  # track the progress:
  counter += 1
  if counter % 1000000 == 0:
    print str(counter), "lines processed"

# write last fragments that are in memory
if GTblock:
  writeBlock(GTblock, fileoutput)
if GTbefore:
  writeBefore(GTbefore, fileoutput)

GTfile.close()
phasedFile.close()
fileoutput.close()

print('Done!')
