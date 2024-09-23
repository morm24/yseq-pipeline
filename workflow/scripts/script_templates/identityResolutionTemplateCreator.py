# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:13:00 2020

@author: hunte
"""

#pip install
import os
import psutil
from multiprocessing.pool import ThreadPool as Pool

import sys
import time
import subprocess
    
def createBED(pos, bedFile):
    with open(bedFile, "w") as w:
        w.write("\t".join(["chrY",str(pos-1-500),str(pos + 500)]) + "\n")
    w.close()


    
def getSequenceFromFasta(bedFile, referenceFile):
    sequence = subprocess.check_output(["bedtools", "getfasta", "-fi", referenceFile, "-bed", bedFile]).decode().split("\n")[1]
    return sequence

def getBLAT(userSeq, blatFile):        
    a = subprocess.check_output(['curl', '-d', "userSeq=" + userSeq + "&type=DNA", '-X', 'POST', 'https://genome.ucsc.edu/cgi-bin/hgBlat']).decode()
    with open(blatFile, "w") as w:
        w.writelines(a)
    w.close

def executeMinimap2(referenceFile, fastQFile, pafOutputFile):
    a = subprocess.check_output(['minimap2', '-c', '-X', referenceFile, fastQFile]).decode()
    with open(pafOutputFile, "w") as w:
        w.writelines(a)
    w.close()

ignoreSequences = ["fix","alt"]

def parsePAF(pos, pafOutputFile):
    fails = {}
    startPosition = pos - 501
    endPosition = pos + 500
    with open(pafOutputFile, "r") as r:
        lines = r.readlines()
        for line in lines:
            splitRow = line.split("\t")
            querySeqLength = int(splitRow[1])
            targSeqName = splitRow[5]
            targStart = int(splitRow[7])
            targEnd = int(splitRow[8])
            residueMatches = int(splitRow[9])
            alignmentBlockLength = int(splitRow[10])
            if any(s in targSeqName for s in ignoreSequences) == False:
                if targStart != startPosition or targEnd != endPosition:
                    minDenom = min([querySeqLength, alignmentBlockLength])
                    percent = float(residueMatches) / float(minDenom)
                    if querySeqLength > 500 and alignmentBlockLength > 500 and percent > 0.95:
                        fails[percent] = str(round(percent * 100,1)) + "% " + targSeqName + ":" + str(targStart) + ".." + str(targEnd) + " " + str(minDenom)

    if len(fails) > 0:
        maxPercent = max(fails.keys())
        return fails[maxPercent]
    return "ok"
        
def parseBLAT(pos, blatOutputFile):
    import re
    startPosition = pos - 500
    fails = []
    with open(blatOutputFile, "r") as r:
        lines = r.readlines()
        foundResults = False
        parseString = "YourSeq"#"chrY:" + str(pos-501) + "-" + str(pos+500)
        for line in lines:
            if foundResults:
                if parseString in line:
                    resultsSplit = line.split(parseString)[2]
                    stripped = re.sub(' +', ' ', resultsSplit).strip().split(" ")
                    targStartPos = int(stripped[7])
                    targEndPos = int(stripped[8])
                    targSeqName = stripped[5]
                    if any(s in targSeqName for s in ignoreSequences) == False:

                        if targStartPos != startPosition or targEndPos != startPosition + 1000:
                            percent = float(stripped[4].replace("%","")) / 100
                            targStart = int(stripped[1])
                            targEnd = int(stripped[2])
                            span = int(stripped[9])
    
                            if targEnd - targStart > 500 and span > 500 and percent >= .95:
                                fails.append(str(round(percent * 100,1)) + "% " + targSeqName + ":" + str(targStart) + ".." + str(targEnd) + " " + str(min([span,targEnd - targStart])))
            else:
                if "----------" in line:
                    foundResults = True
    return fails

def createFastQ(seq, fastQFilename):    
    with open(fastQFilename, "w") as w:
        w.write("@blah\n")
        w.write(seq + "\n")
        w.write("+\n")
        dummy = "".join("F" for f in seq)
        w.write(dummy + "\n")
    w.close()
    
def oneOff(referenceFile, refIndexFile, position):
    bedFile = "bed.bed"+str(position)
    fastQFile = "fastq.fastq"+str(position)
    
    createBED(position, bedFile)
    seq = getSequenceFromFasta(bedFile, referenceFile)
    createFastQ(seq, fastQFile)
    pafOutputFile = "alignment.paf"+str(position)
    executeMinimap2(refIndexFile, fastQFile, pafOutputFile)
    return parsePAF(position, pafOutputFile)

def testPAF(refIndexFile, position):
    pafOutputFile = "alignment.paf"+str(position)
    return parsePAF(position, pafOutputFile)

def blatOneOff(referenceFile, position):
    bedFile = "bed.bed"
    
    createBED(position, bedFile)
    blatFile = str(position) + "_BLAT"
    seq = getSequenceFromFasta(bedFile, referenceFile)
    print(seq)
    getBLAT(seq, blatFile)
    fails = parseBLAT(position, blatFile)
    if len(fails) > 0:
        return fails[0]
    else:
        return "ok"

def checkRanges(pos, refseq):
    rangeMap = {"PAR1": [0,2781478],
                "CEN": [10072349,11686749],
                "DYZ19": [20054913, 20351053],
                "PostPali":[26637970,26673209],
                "PAR2":[26673210,57227415]}
    if refseq == "hs1":
        rangeMap = {"PAR1": [0,2458319],
                    "CEN": [10391466,12593749],
                    "DYZ19": [20961438, 21226267],
                    "PAR2":[62454570,62460029]}
    for r in rangeMap:
        if pos >= rangeMap[r][0] and pos <= rangeMap[r][1]:
            return r
    return "ok"
    

def lineProcessing(line, w, referenceFile, refIndexFile, refseq, output, passing):
    rowSplit = line.replace("\n","").split("\t")
    pos = rowSplit[1]

    #split row split again, so extraction of info field and DP4 values are possible
    info = rowSplit[7]
    infoSplit = info.split(";")
    cov = 0

    #calculate coverage at the snp position
    for dp4 in infoSplit:
        if (dp4.find("DP4") +1):     #if found: value > 0 == true. if not found: -1+1=0 == false
            dp4 = dp4.replace("DP4=","")
            dp4Split = dp4.split(",")
            cov = int(dp4Split[0]) + int(dp4Split[1]) + int(dp4Split[2]) + int(dp4Split[3])
    #print("coverage\n" + str(cov))
    idVal = "chrY:" + pos
    intpos = int(pos)
    rowSplit[2] = idVal
                
    checkedRange = checkRanges(intpos, refseq)
    if checkedRange == "ok":
        if int(cov) == 2:
            checkedRange = "2x"
            passing.append(pos)
    #print(cov)
    if checkedRange == "ok":
        checkedRange = oneOff(referenceFile, refIndexFile, intpos)

        os.remove("fastq.fastq"+str(pos))
        os.remove("bed.bed"+str(pos))
        os.remove("alignment.paf"+str(pos))

    if checkedRange == "ok":
        passing.append(pos)
        

    blatLink = "https://ket.yseq.de:8443/Finch/FTDNA_EDI/Primer3/input?seq=ChrY:" + str(intpos - 500) + ".." + str(intpos+500) + "&name=" + idVal
    uscsGenomeBrowser = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" + refseq + "&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrY%3A" + str(intpos - 5) + "%2D" + str(intpos + 5) + "&hgsid=661900165_RCnZU2Sd4dWMqySzhwCkmOqVSxQI"
    rowSplit = rowSplit + [checkedRange, blatLink, uscsGenomeBrowser]

    
    output.append("\t".join(rowSplit) + "\n")



def analyzeNovelSNPs(vcfFile, outputFile, referenceFile):
    
    #calculate num_threads by max possible memory usage (8Gb usage per thread)
    maxCores = psutil.virtual_memory().total // (1024 ** 3) // 8
    num_threads = os.cpu_count() if maxCores > os.cpu_count() else maxCores

    print("we are using " + str(num_threads) + " threads")

    refIndexFile = referenceFile + ".mmi"
    refseq = "hg38"
    if referenceFile.find("hs1") > 0:
        refseq = "hs1"
    start_time = time.time()
    passing = []

    output = []

    #open read- and write-files (tsv)
    with open(vcfFile, "r") as f:
        with open(outputFile, "w") as w:
            lines = f.readlines()
            headers = lines[0].replace("\n","").split("\t")
            headers = headers + ["Exclusions","BLAT link","UCSC Genome Browser"]
            w.write("\t".join(headers) + "\n")



            #create multi threading pool
            pool = Pool(num_threads)

            for line in lines[1:]:
                #start subfunction multi threaded: #lineProcessing(line,w,referenceFile, refIndexFile, refseq, output, passing)
                pool.apply_async(lineProcessing,(line,w,referenceFile, refIndexFile, refseq, output, passing,))
            pool.close()
            pool.join()   

            #sort unsorted output array by collumn 2 (start pos)
            output = sorted(output, key=lambda x: int(x.split('\t')[1]))
            print("\n\n\n\n\nfirst line of array \n\n\n\n")
            if output: 
                print(output[0])
            else:
                print("empty list")
            print("\n\n\n\n\n")
            for text in output:
                #print(text)
                w.write(text)
        w.close()
    f.close()

    with open("novelPassingPositionsForBLATCheck.txt", "w") as w:
        w.write(",".join(passing))
    w.close()    
    end_time = time.time()
    print("elapsed time: " + str(end_time - start_time) + "s")
          
from time import sleep

if len(sys.argv) > 3:
    mode = sys.argv[1]
    if mode == "-batch":
        vcfFile = sys.argv[2]
        outputFile = sys.argv[3]
        referenceFile = sys.argv[4]
        analyzeNovelSNPs(vcfFile, outputFile, referenceFile)        
    else:
        if mode == "-blat":
            delay = 16
            results = []
            referenceFile = sys.argv[2]
            positions = []
            oktotal = 0
            for pos in sys.argv[3].split(","):
                positions.append(int(pos))
            for pos in positions:
                result = blatOneOff(referenceFile, pos)
                if result == "ok":
                    oktotal = oktotal + 1
                results.append(str(pos) + " " + result)
                sleep(delay)  
            for result in results:
                print(result)
            print(str(oktotal) + " of " + str(len(positions)) + " pass")
        else:
            if mode == "-minimap":
                referenceFile = sys.argv[2]
                position = int(sys.argv[3])
                print(oneOff(referenceFile, position))
