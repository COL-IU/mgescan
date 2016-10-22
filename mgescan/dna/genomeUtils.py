from RBTree import *
import re
import math
import sys

class Gene:
 def __init__(self,start,end,strand,index,seq):
  self.start = start
  self.end = end
  self.strand = strand # 0 for forward and 1 for reverse
  self.m1 = None
  self.m2 = None
  self.m3 = None
  self.info = None
  self.index = index
  self.seq = seq
  
 #def __cmp__(self, other):
 # if other == None:
 #  return False
 # return self.start == other.start and self.end == other.end and self.strand == other.strand
  
 def __hash__(self):
  return hash(str(self.start)+" "+str(self.end)+" "+str(self.strand))
  
 def write(self):
  print str(self.start) + "\t" + str(self.end),
  if self.strand == 1:
   print "\t-\t" + str(self.getFrame()-2),
  else:
   print "\t+\t" + str(self.getFrame()+1),
  print "" 
  
 def toString(self):
  retString = str(self.start) + "\t" + str(self.end)
  if self.strand == 0:
   retString = retString + "\t+\t" + str(self.getFrame()+1) + "\t" + str(self.info)
  else:
   retString = retString + "\t-\t" + str(self.getFrame()+1) + "\t" + str(self.info)
  return retString
	
 def getFrame(self):
  if self.seq == "":
   return -2
  seq = self.seq
  if seq[-3:] == "TAG" or seq[-3:] == "TGA" or seq[-3:] == "TAA":
   if self.strand == 0:
    return self.end%3
   else:
    return ((self.start-1)%3)+3
  elif seq[:3] == "ATG" or seq[:3] == "GTG" or seq[:3] == "TTG":
   if self.strand == 0:
    return (self.start-1)%3
   else:
    return (self.end%3)+3
  else:
   return -1
   
def compareFrames(frame1,frame2):
 if frame1 == -1 or frame2 == -1:
  return False
 return frame1 == frame2

def testPrint():
 print "hi"
 
def proteinToDNAInterval(pep_start,pep_end,DNALength,frame):
 if frame > 3:
  frame = frame - 3
  nuc_start = ((pep_start-1)*3)+((DNALength-frame+1)%3)
  nuc_end = ((pep_end-1)*3)+((DNALength-frame+1)%3)
  return [max(0,(DNALength-nuc_end)-3)+1,(DNALength-nuc_start),-1]
 else:
  return [((pep_start-1)*3)+frame,((pep_end-1)*3)+frame+2,1]

#~/softwares/cross_match <input1> <input2> -tags 1> outToParse 2> /dev/null
def readCrossMatch(filename):
 f = open(filename)
 out = []
 for line in f:
  words = line.split(None)
  if len(words) == 0:
   continue
  if words[0] == "ALIGNMENT":
   out.append(words)
 f.close()
 return out

def readGFF3(filename,seqFile):
 if seqFile != None:
  seqDict = readFASTAinDict(seqFile)
 else:
  seqDict = {}
 fastalines = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 for line in lines:
  if line[0] == "#":
   continue
  words = line.split(None)
  header = words[0]
  start = int(words[3])
  end = int(words[4])
  if header in seqDict:
   seq = seqDict[header][start-1:end]
  else:
   seq = ""
  if words[6] == '+':
   strand = 0
  else:
   strand = 1
   seq = reverse(complement(seq))
  gene = Gene(start,end,strand,0,seq)
  gene.info = line.strip()
  if header in fastalines.keys():
   fastalines[header][0].append([start,end,gene])
   fastalines[header][1].append(gene)
  else:
   fastalines[header] = [[[start,end,gene]],[gene]]
 for header in fastalines.keys():
  fastalines[header][0] = RBTree(fastalines[header][0])
 return fastalines

def readFragOutFmt(filename,seqFile):
 if seqFile != None:
  seqDict = readFASTAinDict(seqFile)
 else:
  seqDict = {}
 fastalines = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 geneNodes = []
 geneObjects = []
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    line = line.rstrip()
    if header != "":
     fastalines[header] = [RBTree(geneNodes),geneObjects]
     geneNodes = []
     geneObjects = []
    header = line.strip()[1:]
   else:
    words = line.split(None)
    start = int(words[0])
    end = int(words[1])
    if header in seqDict:
     seq = seqDict[header][start-1:end]
    else:
     seq = ""
    if words[2] == '+' or words[2] == '.':
     strand = 0
    else:
     strand = 1
     seq = reverse(complement(seq))
    gene = Gene(start,end,strand,len(geneNodes),seq)
    gene.info = line.strip()
    geneNodes.append([start,end,gene])
    geneObjects.append(gene)
 fastalines[header] = [RBTree(geneNodes),geneObjects]
 return fastalines

def complement(sequence):
#This function defines the nucleotide complement dictionary and utilizes it to return a complementary DNA sequence.
 dict = {"A":"T", "T":"A", "C":"G", "G":"C", "X":"X", "O":"O", "?":"?", "N":"N", "a":"t", "t":"a", "c":"g", "g":"c"}
 complement = ""
 for i in range(0,len(sequence)):
  if sequence[i] in dict.keys():
   complement = complement + dict[sequence[i]]
  else:
   complement = complement + sequence[i]
 return complement

def reverse(sequence):
#This function simply returns the reverse of a sequence. This function coupled with the above complement function
#will give a reverse complementary DNA sequence.
 return sequence[::-1]

def readRegulonFile(regFile,seqFile):
 singleSeqDs = readFASTA(seqFile)
 genome = singleSeqDs[0][1]
 geneNodes = []
 geneObjects = []
 f = open(regFile)
 lines = f.readlines()
 f.close()
 for i in range(1,len(lines)):
  line = lines[i]
  words = line.split(",")
  if len(words) > 2:
   start = int(words[2])
   end = int(words[3])
   if words[6] == 'forward':
    strand = 0
    seq = genome[start-1:end]
   else:
    strand = 1
    seq = reverse(complement(genome[start-1:end]))
   gene = Gene(start,end,strand,len(geneNodes),seq)
   gene.info = line.strip()
   geneNodes.append([start,end,gene])
   geneObjects.append(gene)
 return [RBTree(geneNodes),geneObjects] 
 
def readPttFile(pttFile,seqFile):
 singleSeqDs = readFASTA(seqFile)
 genome = singleSeqDs[0][1]
 geneNodes = []
 geneObjects = []
 f = open(pttFile)
 lines = f.readlines()
 f.close()
 for i in range(3,len(lines)):
  line = lines[i]
  words = line.split("\t")
  if len(words) > 2:
   startEnd = words[0].split("..")
   start = int(startEnd[0])
   end = int(startEnd[1])
   if words[1] == '+':
    strand = 0
    seq = genome[start-1:end]
   else:
    strand = 1
    seq = reverse(complement(genome[start-1:end]))
   gene = Gene(start,end,strand,len(geneNodes),seq)
   gene.info = line.strip()
   geneNodes.append([start,end,gene])
   geneObjects.append(gene)
 return [RBTree(geneNodes),geneObjects]

def getAllMatches(gene,tree,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()
 matchNodes = tree.intervalSearchAll([queryStart,queryEnd])
 matchedGenes = []
 for matchNode in matchNodes:
  matchedGenes.append(matchNode.obj)
 matchedGenes.sort(key=lambda x: x.start, reverse=False)
 return matchedGenes

def getAnyMatch(gene,tree,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()
 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 else:
  return matchedGene

def getMatch(gene,tree,overlapThreshold,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()

 overlapPercent = 0.0
 if matchedStart != 0 and matchedEnd != 0:
  overlapNumber = 0.0 + min(matchedEnd,queryEnd)-max(matchedStart,queryStart)
  queryPercent = (overlapNumber/(queryEnd-queryStart))*100
  matchPercent = (overlapNumber/(matchedEnd-matchedStart))*100
  if queryPercent > overlapThreshold and matchPercent > overlapThreshold and (compareFrames(queryFrame,matchedFrame) or ignoreFrame):
   return [matchedGene,queryPercent,matchPercent]
 return -1

def getMatchQuery(gene,tree,overlapThreshold,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()

 overlapPercent = 0.0
 if matchedStart != 0 and matchedEnd != 0:
  overlapNumber = 0.0 + min(matchedEnd,queryEnd)-max(matchedStart,queryStart)
  overlapPercent = (overlapNumber/(queryEnd-queryStart))*100
  if overlapPercent > overlapThreshold and (compareFrames(queryFrame,matchedFrame) or ignoreFrame):
   return [matchedGene,overlapPercent]
 return -1

def getMatchSubject(gene,tree,overlapThreshold,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()

 overlapPercent = 0.0
 if matchedStart != 0 and matchedEnd != 0:
  overlapNumber = 0.0 + min(matchedEnd,queryEnd)-max(matchedStart,queryStart)
  overlapPercent = (overlapNumber/(matchedEnd-matchedStart))*100
  if overlapPercent > overlapThreshold and (compareFrames(queryFrame,matchedFrame) or ignoreFrame):
   return [matchedGene,overlapPercent]
 return -1

'''
def getLargerMatch(gene,tree,overlapThreshold):
 overlap = True
 start = gene.start
 end = gene.end
 frame = gene.frame

 match = tree.intervalSearch([start,end])

 matchedGene = match.key
 low = matchedGene["low"]
 high = matchedGene["high"]
 matchFrame = match.frame

 if frame == matchFrame:
  overlapPercent = 0.0
  if low != 0 and high != 0:
   overlapNumber = 0.0 + min(high,end)-max(low,start)
   overlapPercent = (overlapNumber/min((high-low),(end-start)))*100
   if overlapPercent > overlapThreshold:
    return match.info[0]
 return -1
'''
def getExactMatch(gene,tree,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()
 if matchedStart != 0 and matchedEnd != 0:
  if matchedStart == queryStart and matchedEnd == queryEnd and (compareFrames(queryFrame,matchedFrame) or ignoreFrame):
   return matchedGene
 return -1

def readClustalW(filename):
 data = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 header = ""
 for i in range(1,len(lines)):
  line = lines[i]
  words = line.split(None)
  if len(words) != 0:
   if words[0] not in data.keys():
    if header == "":
     header = words[0]
    data[words[0]] = words[1].strip()
   else:
    data[words[0]] = data[words[0]] + words[1].strip()
 return [header,data] 

def readFASTA(filename,upper=False):
 fastalines = []
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 seq = ""
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    if seq != "":
     fastalines.append([header,seq])
     seq = ""
    header = line.strip()[1:]
   else:
    line = re.sub(r'\s','',line)
    if upper:
     line = line.upper()
    seq = seq + line
 if seq != "":
  fastalines.append([header,seq])
 return fastalines

def readFASTAinDict(filename,upper=False):
 fastalines = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 seq = ""
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    if seq != "":
     words = header.split(None)
     fastalines[words[0][1:]] = seq
     seq = ""
    header = line
   else:
    line = re.sub(r'\s','',line)
    if upper:
     line = line.upper()
    seq = seq + line
 if seq != "":
  words = header.split(None)
  fastalines[words[0][1:]] = seq
 return fastalines

def readEMBLinDict(filename):
 out = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0:2] != "ID":
  sys.exit("Invalid EMBL file: "+filename)
 data = []
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0:2]=="ID":
    if len(data) > 0:
     words = header.split(None)
     out[words[1]] = data
     data = []
    header = line
   data.append(line)
    
 if len(data) > 0:
  words = header.split(None)
  out[words[1]] = data
 return out

#Define a function to print the matrices in readable form (For debugging purposes)
def printMatrix(matrix,seq1,seq2):
	columns = len(seq1)+1
	rows = len(seq2)+1
	print "#\t*\t",
	for i in range(0,columns-1):
		print seq1[i]+"\t",
	print "\n"
	for i in range(0,rows):
		if i == 0:
			print "*\t",
		else:
			print seq2[i-1]+"\t",
		for j in range(0,columns):
			print str(matrix[i][j])+"\t",
		print "\n"

#Define a function to do the Smith-Waterman alignment of two sequences
def SWalign(matrix,seq1,seq2,gap_init,gap_ext,match,mismatch):
 n = len(seq1) #Length of first sequence
 m = len(seq2) #Length of second sequence
 sigma = gap_ext #Gap extension penalty
 rho = gap_init #Gap initiation penalty
 S = [ [0]*(n+1) for i in range(m+1) ] #Score matrix
 I = [ [0]*(n+1) for i in range(m+1) ] #Insertion score matrix
 D = [ [0]*(n+1) for i in range(m+1) ] #Deletion score matrix
 T = [ [0]*(n+1) for i in range(m+1) ] #Traceback matrix
 max_score = 0
 max_i = 0
 max_j = 0
 if n != 0 and m != 0:
  #Initialize the matrices (first row/column)
  for i in range(0,n+1):
   D[0][i] = 0
   I[0][i] = -(11+i)
   S[0][i] = 0
  for i in range(0,m+1):
   D[i][0] = -(11+i)
   I[i][0] = 0
   S[i][0] = 0
  #Fill up the rest of the cells by calculating the scores
  for i in range(1,m+1):
   for j in range(1,n+1):
    D[i][j] = max(D[i-1][j]+sigma, S[i-1][j]+sigma+rho)
    I[i][j] = max(I[i][j-1]+sigma, S[i][j-1]+sigma+rho)
    if matrix == None:
     if seq1[j-1] == seq2[i-1]:
      a = S[i-1][j-1] + match
     else:
      a = S[i-1][j-1] + mismatch
    else:
     a = S[i-1][j-1]+matrix[seq1[j-1]+seq2[i-1]]
    b = I[i][j]
    c = D[i][j]
    S[i][j] = max(0, a, b, c)
    
    #Markers for tracing back
    if S[i][j] == a:
     T[i][j] = 1 #Match/Mismatch
    elif S[i][j] == b:
     T[i][j] = 2 #Insertion
    elif S[i][j] == c:
     T[i][j] = 3 #Deletion
    
    #Keep track of the maximum score and its indices
    if S[i][j] > max_score:
     max_score = S[i][j]
     max_i = i
     max_j = j

  #printMatrix(D,seq1,seq2)
  #printMatrix(I,seq1,seq2)
  #printMatrix(S,seq1,seq2)	  
  #print max_score
  score = max_score	 
  #Traceback
  s1 = ""
  s2 = ""
  e1 = max_j
  e2 = max_i
  while max_score != 0:
   #Match/Mismatch
   if T[max_i][max_j] == 1:
    s1 = s1 + seq1[max_j-1]
    s2 = s2 + seq2[max_i-1]
    max_i = max_i-1
    max_j = max_j-1
   #Deletion
   elif T[max_i][max_j] == 3:
    s1 = s1 + "-"
    s2 = s2 + seq2[max_i-1]
    max_i = max_i-1
   #Insertion
   elif T[max_i][max_j] == 2:
    s1 = s1 + seq1[max_j-1]
    s2 = s2 + "-"
    max_j = max_j-1
   
   max_score = S[max_i][max_j]

  b1 = max_j+1
  b2 = max_i+1
  s1 = s1[::-1]
  s2 = s2[::-1]
  return [s1,s2,b1,e1,b2,e2,score]
 else:
  return None
