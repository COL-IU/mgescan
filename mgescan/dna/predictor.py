import sys
from Bio.Blast import NCBIXML
from genomeUtils import *
from scriptUtils import *

def top2(lst):
 slst = sorted(lst,key=lambda x: x[0], reverse=True)
 result = []
 top2Keys = set()
 for item in slst:
  if len(top2Keys) < 2 or item[0] in top2Keys:
   result.append(item)
   top2Keys.add(item[0])
  else:
   break
 return result

def getBlastAlignments(seq1,seq2,tempdir):
 temp = open(tempdir+"query","w")
 temp.write(">query\n"+seq1+"\n")
 temp.close()
 temp = open(tempdir+"subject","w")
 temp.write(">subject\n"+seq2+"\n")
 temp.close()
 cmd = " ".join(['blastn','-query',tempdir+"query",'-subject',tempdir+"subject",'-out',tempdir+"results.xml",'-outfmt','5','-dust','no','-word_size','5','-strand','plus','-evalue','0.15','-reward','1','-penalty','-3','-gapopen','5','-gapextend','2'])
 runCommand(cmd,None)
 alignments = []
 if os.path.getsize(tempdir+"results.xml") == 0:
  return []
 result_handle = open(tempdir+"results.xml")
 blast_records = NCBIXML.parse(result_handle)
 blast_records = list(blast_records)
 for blast_record in blast_records:
  for alignment in blast_record.alignments:
   for hsp in alignment.hsps:
    if hsp.frame[0] >= 0 and hsp.frame[1] >= 0:
     alignments.append([hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.score])
 return alignments

def getBestBlastScore(seq,db,tempdir):
 temp = open(tempdir+"tir","w")
 temp.write(">query\n"+seq+"\n")
 temp.close()
 cmd = " ".join(['blastn','-query',tempdir+"tir",'-db',db,'-out',tempdir+"tirResults.tbl",'-outfmt','\"6 score\"','-dust','no','-word_size','4','-max_target_seqs','1','-reward','3','-penalty','-2','-gapopen','5','-gapextend','5'])
 runCommand(cmd,None)
 if os.path.getsize(tempdir+"tirResults.tbl") == 0:
  return 0
 result_handle = open(tempdir+"tirResults.tbl")
 return float(result_handle.readline())

def loadModels(paramdir):
 models = {}
 for filename in listFiles(paramdir, '*.param'):
  name = removeExt(filename)
  loading = {}
  f = open(paramdir+filename) 
  for line in f:
   words = line.split(":")
   loading[words[0]] = float(words[1].strip())
  models[name] = loading
  f.close()
 return models

def scoreCandidate(seq,model,candidate):
 nucs = ['A','C','G','T']
 a = candidate[0]
 b = candidate[1]
 c = candidate[2]
 d = candidate[3]
 e = candidate[4]
 f = candidate[5]
 lp1 = 0.0
 counter = 0
 for i in range(a-1,f-1):
  lp = 0
  if seq[i] in nucs and seq[i+1] in nucs:
   if (i >= b and i < c) or (i >= d and i < e):
    if "0->"+seq[i:i+2] in model.keys():
     lp = math.log(model["0->"+seq[i:i+2]],2)
     counter += 1
    else:
     sys.exit("The model parameters file is invalid")
   elif i < b or i >= e:
    if "1->"+seq[i:i+2] in model.keys():
     lp = math.log(model["1->"+seq[i:i+2]],2)
     counter += 1
    else:
     sys.exit("The model parameters file is invalid")
   else:
    if "2->"+seq[i:i+2] in model.keys():
     lp = math.log(model["2->"+seq[i:i+2]],2)
     counter += 1
    else:
     sys.exit("The model parameters file is invalid")
  else:
   lp = math.log(0.0625,2)
   counter += 1
  lp1 += lp

 lp0 = counter * math.log(0.0625,2)
 score = lp1 - lp0
 return score
 
def runHmmSearch(fastaPath,hmmPath,outDir):
 cmd = " ".join(['hmmsearch','-E','0.0001','--noali','--nobias','--domtblout',outDir+"domtbl",hmmPath,fastaPath])
 runCommand(cmd)
 results = []
 f = open(outDir+"domtbl")
 for line in f:
  if len(line) > 0:
   if line[0] == '#':
    continue
   words = line.split(None)
   seqName = words[0]
   hmmName = words[3]
   qlen = int(words[5])
   ievalue = float(words[12])
   score = float(words[13])
   bias = float(words[14])
   start = int(words[19])
   end = int(words[20])
   acc = float(words[21])
   alilen = int(words[20]) - int(words[19]) + 1
   aliPercent = (alilen+0.0)/qlen
   results.append([seqName,hmmName,qlen,score,bias,start,end,aliPercent,ievalue,acc])
 f.close()
 myLog("Finished hmmsearch.")
 return results
 
def readHmmerResults(domtblPath):
 myLog("Reading HMMER results from "+domtblPath)
 results = []
 f = open(domtblPath)
 for line in f:
  if len(line) > 0:
   if line[0] == '#':
    continue
   words = line.split(None)
   seqName = words[0]
   hmmName = words[3]
   qlen = int(words[5])
   ievalue = float(words[12])
   score = float(words[13])
   bias = float(words[14])
   start = int(words[19])
   end = int(words[20])
   acc = float(words[21])
   alilen = int(words[20]) - int(words[19]) + 1
   aliPercent = (alilen+0.0)/qlen
   results.append([seqName,hmmName,qlen,score,bias,start,end,aliPercent,ievalue,acc])
 f.close()
 return results
 
def translateSeqs(DNAPath,proteinPath):
 cmd = " ".join(['transeq','-sequence',DNAPath,'-outseq',proteinPath,'-frame=6'])
 runCommand(cmd)

#################################################################################################
#################################################################################################
#################################################################################################
 

def runDNA(indir,outdir,basepath):

 myLog("Beginning DNA transposon prediction...")

 indir = cleanDirPath(indir)
 outdir = cleanDirPath(outdir)
 dnaoutdir = cleanDirPath(outdir+"dna/")
 makeDir(dnaoutdir)
 paramdir = basepath + "/dna/params/"
 tirdir = basepath + "/dna/TIRDB/"

 tempdir = makeTempDir(dnaoutdir)
 
 cmd = "cat "+indir+"* > "+tempdir+"allSeqs.fa"
 runCommandShell(cmd,None)
 
 infile = tempdir+"allSeqs.fa"

 slimits = {"academ":6000,"cmc":15000,"ginger":6000,"harbinger":9000,"hat":6000,"kolobok":6000,"mariner":6000,"merlin":3000,"mule":15000,"novosib":6000,"p":6000,"piggybac":6000,"sola1":3000,"sola2":3000,"sola3":5000,"transib":3000,"zator":3000,"1360":6000,"hopper":1000}

 hmmpath = basepath+"/dna/domains.hmm"
 pepfile = tempdir+"input.pep"

 seqDict = readFASTAinDict(infile,True)
 models = loadModels(paramdir)

 lengths = {}
 for key in seqDict.keys():
  lengths[key] = len(seqDict[key])

 if os.path.isfile(dnaoutdir+"domtbl") == False:
  translateSeqs(infile,pepfile)
  hmmerResults = runHmmSearch(pepfile,hmmpath,dnaoutdir)
 else:
  hmmerResults = readHmmerResults(dnaoutdir+"domtbl")

 anchorsDict = {}
 for result in hmmerResults:
  [seqName,hmmName,qlen,score,bias,start,end,aliPercent,ievalue,acc] = result
  if score > 1 and ievalue < 0.01:
   header = seqName[:-2]
   frame = int(seqName[-1])
   length = lengths[header]
   [nuc_start,nuc_end,strand] = proteinToDNAInterval(start,end,length,frame)
   if header in anchorsDict.keys():
    anchorsDict[header].append([hmmName,nuc_start,nuc_end,strand,score,bias,aliPercent,ievalue,acc])
   else:
    anchorsDict[header] = [[hmmName,nuc_start,nuc_end,strand,score,bias,aliPercent,ievalue,acc]]

 dtPred = {}
 for header in seqDict.keys():
  myLog("Analyzing transposase domains in " + header)
  if header not in anchorsDict.keys():
   dtPred[header] = []
  else:
   anchors = anchorsDict[header]
   seq = seqDict[header]
   length = lengths[header]
   groups = []
   num = len(anchors)
   counter = 0
   for [hmmName,nuc_start,nuc_end,hmmStrand,hmmScore,hmmBias,hmmPercent,ievalue,acc] in anchors:
    counter += 1

    #sys.stdout.write(str(counter)+" of "+str(num) + "\t" + hmmName + "\n")

    lstart = max(0,nuc_start-slimits[hmmName])
    lend = nuc_start
    rstart = nuc_end
    rend = min(length,nuc_end+slimits[hmmName])
    seq1 = seq[lstart:lend]
    seq2 = reverse(complement(seq[rstart:rend]))
    group = []
    if seq1 != "" and seq2 != "":
     alignments = getBlastAlignments(seq1,seq2,tempdir)
     for i in range(0,min(3,len(alignments))):
      alignment = alignments[i]
      a = lstart+alignment[0]
      b = lstart+alignment[1]
      c = nuc_start
      d = nuc_end
      e = rend-alignment[3]
      f = rend-alignment[2]
      candidate = [a,b,c,d,e,f]
      myScore = scoreCandidate(seq,models[hmmName],candidate)
      tirScore1 = alignment[4]
      predLength = f-a+1
      tirScore2 = getBestBlastScore(seq[a-1:b],tirdir+hmmName,tempdir)
      totalScore = hmmScore + tirScore1 + tirScore2
      if myScore > 0:
       group.append([totalScore,a,b,c,d,e,f,myScore,hmmScore,tirScore1,tirScore2,hmmName,hmmPercent,ievalue,acc])
    if len(group) > 0:
     topTwo = top2(group)
     groups.append([topTwo[0][1],topTwo])
   
   if len(groups) == 0:
    dtPred[header] = []
   else:
    sortedGroups = sorted(groups,key=lambda x: x[0]) 
    filtered = []

    prev = sortedGroups[0][1]
    filtered.append(prev)
    for i in range(1,len(sortedGroups)):
     group = sortedGroups[i][1]
     if group[0][1] <= prev[0][6]: #overlapping
      if group[0][0] > prev[0][0]:
       filtered.pop()
       filtered.append(group)
       prev = group
     else:
      filtered.append(group)
      prev = group
    dtPred[header] = filtered

 predGenesDict = {}
 o1 = file(dnaoutdir+"original.out","w")
 o2 = file(dnaoutdir+"extended.out","w")
 for header in seqDict.keys():
  o1.write(">"+header+"\n")
  o2.write(">"+header+"\n")
  predGenesList = []
  groups = dtPred[header]
  for group in groups:
   for i in range(0,len(group)):
    [totalScore,a,b,c,d,e,f,myScore,hmmScore,tirScore1,tirScore2,hmmName,hmmPercent,ievalue,acc] = group[i]
    #tirLen = b-a+1
    outString = myljust(a,12) + myljust(f,12) + "+  " + myljust(hmmName,11) + myljust(1,6) + myljust(b-a+1,8) + myljust(c-a+1,8) + myljust(d-a+1,8) + myljust(e-a+1,8) + myljust(f-a+1,8) + myljust(hmmScore,8) + myljust(tirScore1,8) + myljust(tirScore2,8) + "\n"
    if i == 0:
     o1.write(outString)
     gene = Gene(int(a),int(f),0,0,"")
     gene.info = outString
     predGenesList.append(gene)
    o2.write(outString)
  predGenesDict[header] = predGenesList

 o1.close()
 o2.close()

 filteredGenesDict = {}
 ###### Subtract LTRs ######
 ltrs = readGFF3(outdir+"ltr/ltr.gff3",None)
 for header in predGenesDict.keys():
  predGenesList = predGenesDict[header]
  filteredGenesList = []
  
  if header not in ltrs.keys():
   ltrTree = RBTree([])
   ltrGenes = []
  else:
   ltrTree = ltrs[header][0]
   ltrGenes = ltrs[header][1]

  for predGene in predGenesList:
   matches = getAllMatches(predGene,ltrTree,True)
   if len(matches) == 0:
    filteredGenesList.append(predGene)
  filteredGenesDict[header] = filteredGenesList
  
 o3 = file(dnaoutdir+"final.out","w")
 for header in filteredGenesDict.keys():
  o3.write(">"+header+"\n")
  filteredGenesList = filteredGenesDict[header]
  for filteredGene in filteredGenesList:
   o3.write(filteredGene.info)
 o3.close()
 
 o4 = file(outdir+"mge.gff3","w")
 ###### Process GFFs ######
 nonltrs = readGFF3(outdir+"info/nonltr.gff3",None)
 for header in seqDict.keys():
  if header in ltrs.keys():
   ltrList = ltrs[header][1]
  else:
   ltrList = []
  if header in nonltrs.keys():
   nonltrList = nonltrs[header][1]
  else:
   nonltrList = []
  ltrList.extend(nonltrList)
  dnaList = filteredGenesDict[header]
  for dnaItem in dnaList:
   words = dnaItem.info.split(None)
   id = header + "_" + words[0]
   name = words[3]
   dnaItem.info = "\t".join([header, "MGEScan_DNA", "mobile_genetic_element", str(dnaItem.start), str(dnaItem.end), ".", ".", ".", "ID=" + id + ";name="+name])
   ltrList.append(dnaItem)
  
  sortedList = sorted(ltrList,key=lambda x: x.start)
  for item in sortedList:
   o4.write(item.info + "\n")
   
 o4.close()
   
 
 removeDir(tempdir)

 myLog("Completed")
