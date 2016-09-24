#!/usr/bin/env python2
#
# The script merged phased heterozygous sites with homozygous sites that were skipped during phasing by HapCUT.
# 
# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
#
# python mergePHASEDsnpsGT.py -p {phased-data} -g {GT.table} -o {output} -Np 0.15
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
parser.add_argument('-g', '--gt_table', help = 'name of the GT table file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-Np', '--probability_of_N', help = 'probability of N in homozygotes', type=float, required=True)
args = parser.parse_args()

prN = args.probability_of_N

def hetToNsBlock(g):
  if g[0] != g[1] or g == ['.', '.']:
    g = ['N', 'N']
  return g

def hetToNsBefore(g):
  if g[0] != g[1] or g == ['.', '.']:
    g = ['N', 'N']
  else:
    gtN = numpy.random.binomial(1, prN, 1)
    g = [g, ['N', 'N']][int(gtN)]
  return g

def writeBlock(listToWrire, fileout): 
  if all(x[1][0] == x[1][1] for x in listToWrire):
    for i in listToWrire:
      i[1] = ['N', 'N']
  for el in range(len(listToWrire)):
    listToWrireP = '\t'.join(str(i) for sublist in listToWrire[el] for i in sublist)
    fileout.write("%s\n" % listToWrireP)

def writeBefore(listToWrire, fileout):
  for el in range(len(listToWrire)):
    if listToWrire[el][1][0] != listToWrire[el][1][1] or listToWrire[el][1] == ['.', '.']:
      listToWrire[el][1] = ['N', 'N']
    listToWrireP = '\t'.join(str(i) for sublist in listToWrire[el] for i in sublist)
    fileout.write("%s\n" % listToWrireP)
    
def appendGT(coord, g, block):
  if g[0] != g[1] or g == ['.', '.']:
    gt = ['N', 'N']
  else:
    gt = g
  block.append([coord, g])
  
def readPhasedLine(fileName, stat):
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
  stat = 'block'
  GTblo.append([coord, pGT])
  pL = readPhasedLine(pFile, stat)
  pChr  = pL[0]
  pPos = pL[1]
  pGT = pL[2]
  stat = pL[3]
  return (GTblo, pChr, pPos, pGT, stat)

def writeFragment(stat1, stat, GTbe, GTblo, fileout):
  if stat1 == 'before' and stat == 'block':
    writeBefore(GTbe, fileout)
    GTbe = []
    stat1 = stat
  elif stat == 'before' and stat1 == 'block':
    writeBlock(GTblo, fileout)
    GTblo = []
    stat1 = stat
  return (stat1, GTbe, GTblo)

phasedFile = open(args.phased, 'r')
phasedHeader = phasedFile.readline()
phasedHeaderW = phasedHeader.split()
name = phasedHeaderW[2].split('_')[0]

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
indnumber = GTheader.index(name + '.GT')


counter = 0
for line in GTfile:
  words = line.split()
  GTchr = words[0].split('_')[1]
  GTpos = words[1]
  GTcoord = words[0:2]
  GTgt = words[indnumber].split('/')
  stopChr1 = stopChr
  stopPos1 = stopPos
  if stopChr == 'NA':
    writeBlock(GTblock, fileoutput)
    GTblock = []
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos) and status == 'before':
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) < int(stopPos) and status == 'block':
    GTgtN = hetToNsBlock(GTgt)
    appendGT(GTcoord, GTgtN, GTblock)
  elif (GTchr == stopChr and GTpos == stopPos):
    eqLine = appendPhasedGT(GTblock, GTcoord, phasedGT, phasedFile)
    GTblock = eqLine[0]
    stopChr = eqLine[1]
    stopPos = eqLine[2]
    phasedGT = eqLine[3]
    status = eqLine[4]
    wF = writeFragment(status1, status, GTbefore, GTblock, fileoutput)
    status1 = wF[0]
    GTbefore = wF[1]
    GTblock = wF[2]
  elif int(GTchr) != int(stopChr):
    status = 'before'
    writeBlock(GTblock, fileoutput)
    GTblock = []
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
  elif int(GTchr) == int(stopChr) and int(GTpos) > int(stopPos):
    GTgtN = hetToNsBefore(GTgt)
    appendGT(GTcoord, GTgtN, GTbefore)
    while int(GTchr) >= int(stopChr) and int(GTpos) > int(stopPos):
      phasedL = readPhasedLine(phasedFile, status)
      stopChr  = phasedL[0]
      stopPos = phasedL[1]
      phasedGT = phasedL[2]
    wF = writeFragment(status1, status, GTbefore, GTblock, fileoutput)
    status1 = wF[0]
    GTbefore = wF[1]
    GTblock = wF[2]
  else:
    print('ERROR')
    
  counter += 1
  if counter % 1000000 == 0:
    print str(counter), "lines processed"
  else:
    continue

writeBlock(GTblock, fileoutput)
writeBefore(GTbefore, fileoutput)

GTfile.close()
phasedFile.close()
fileoutput.close()
print('Done!')
