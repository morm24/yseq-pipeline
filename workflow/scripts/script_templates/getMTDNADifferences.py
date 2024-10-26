# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 17:50:06 2020

@author: hunter, moritz
"""

import sys
import subprocess



#name output file
haplogrepOutputFile = "haplogrep.tsv"



def process():
    subprocess.check_output(["java", "-jar", haplogrepJar, "--in", vcfFile, "--format", "vcf", "--extend-report", "--out", haplogrepOutputFile])
    hg = None
    snpsString = None
    with open(haplogrepOutputFile, "r") as r:
        headers = r.readline()
        info = r.readline().split("\t")
        if len(info) > 2:
            hg = info[2].replace("\"","")
            snpsString = info[9].replace("\n","").replace("\"","")
    r.close()
    
    toIgnore = ["3106d"]
    if hg != None and snpsString != None:
        snps = snpsString.split(" ")
        for snp in toIgnore:
            if snp in snps:
                snps.remove(snp)
    else:
        snps = []
    with open(outFile, "w") as w:
        w.write("\t".join(["position","mutation"]) + "\n")
        for snp in snps:
            w.write("\t".join([snp[:-1],snp[-1:]]) + "\n")
    
    w.close()
    print(hg)

def getSNPsFromOutFile(outFile):
    snps = []
    with open(outFile, "r") as r:
        lines = r.readlines()
        for line in lines[1:]:
            tabs = line.replace("\n","").split("\t")
            snps.append(tabs[0] + tabs[1])
    return snps

def getTextFromPhyloEq(phyloEqFile):
    textLines = []
    with open(phyloEqFile, "r") as r:
        lines = r.readlines()
        for line in lines[1:]:
            tabs = line.replace("\n","").split("\t")
            marker = tabs[0]
            mutation = tabs[1]
            call = tabs[2]
            if call == "no call":
                call = " (no call)"
            if call == "-":
                call = "- (!)"
            if mutation != "unknown":
                textLines.append(marker + " " + mutation + call)            
                    
    return "\n".join(textLines)

def getTextFromDownstr(downstrFile):
    textLines = []
    with open(downstrFile, "r") as r:
        lines = r.readlines()
        lastBranch = False
        for line in lines[1:]:
            tabs = line.replace("\n","").split("\t")
            branch = tabs[0]
            marker = tabs[1]
            mutation = tabs[2]
            call = tabs[3]
            if call == "no call":
                call = " (no call)"
            if mutation != "unknown":
                if lastBranch != branch:
                    textLines.append("")
                    textLines.append(branch + ":")
                textLines.append("  " + marker + " " + mutation + call)
                lastBranch = branch
    if len(textLines) == 0:
        textLines = ["None"]
    return "\n".join(textLines)
        
def update(resultsSummaryFileIn, resultsSummaryFileOut, mtDNAFile, phyloEqFile, downstrFile, mtdnaTagStr, phyloEqsTagStr, downstrTagStr, novelSNPsTagStr):
    snps = getSNPsFromOutFile(mtDNAFile)
    mtDNAsnpsString = " ".join(snps)
    
    phyloEqText = getTextFromPhyloEq(phyloEqFile)
    downstrText = getTextFromDownstr(downstrFile)
    
    with open(resultsSummaryFileIn, "r") as r:
        with open(resultsSummaryFileOut, "w") as w:
            lines = r.readlines()
            replaceLine = False
            ignore = False
            for line in lines:
                if replaceLine == "mtDNA":
                    w.write(mtDNAsnpsString + "\n")        
                    replaceLine = False                    
                else:
                    if replaceLine == "phyloeq":
                        w.write("\n" + phyloEqText + "\n")
                        replaceLine = False
                        ignore = True
                    else:
                        if replaceLine == "downstr":
                            w.write(downstrText + "\n")
                            replaceLine = False
                            ignore = True
                        else:                            
                            if mtdnaTagStr in line:
                                replaceLine = "mtDNA"
                                w.write(line)
                            else:
                                if phyloEqTagStr in line:
                                    replaceLine = "phyloeq"
                                    w.write(line)
                                else:
                                    if downstrTagStr in line:
                                        replaceLine = "downstr"
                                        w.write("\n")
                                        w.write(line)
                                    else:
                                        if novelSNPsTagStr in line:
                                            ignore = False
                                            w.write("\n\n")
                                            w.write(line)
                                        else:            
                                            if ignore == False:
                                                w.write(line)
        w.close()
    r.close()

def getSNPsFromPhyloEqNew(phyloEqFile):
    markerAlleles = []
    with open(phyloEqFile, "r") as r:
        lines = r.readlines()
        for line in lines[1:]:
            tabs = line.replace("\n","").split("\t")
            markers = tabs[0].split("/")
            mutation = tabs[1]
            call = tabs[2]
            if mutation != "unknown" and call != "no call":
                for marker in markers:
                    markerAlleles.append((marker, mutation + call))            
    return markerAlleles

def getSNPsFromDownstrNew(downstrFile):
    markerAlleles = []
    with open(downstrFile, "r") as r:
        lines = r.readlines()
        for line in lines[1:]:
            tabs = line.replace("\n","").split("\t")
            markers = tabs[1].split("/")
            mutation = tabs[2]
            call = tabs[3]
            if mutation != "unknown" and call != "no call":
                for marker in markers:
                    markerAlleles.append((marker, mutation + call))            
    return markerAlleles

def addAlleles(mtdnaOutFile, addAllelesFile, sampleName, mtdnaSourceFile, phyloEqFile, downstrFile, ydnaSourceFile):
    snps = getSNPsFromOutFile(mtdnaOutFile)
    print(addAllelesFile)
    with open(addAllelesFile, "w") as w:
        w.write("\t".join(["#Sample Name","Marker Name","Allele","Comments"]) + "\n")
        for snp in snps:
            w.write("\t".join([sampleName,"mtDNA:" + snp[:-1], snp[-1:] + "+", mtdnaSourceFile]) + "\n")
        phyloEqSNPs = getSNPsFromPhyloEqNew(phyloEqFile)
        for markerAllele in phyloEqSNPs:
            w.write("\t".join([sampleName, markerAllele[0], markerAllele[1], ydnaSourceFile]) + "\n")
        downstrSNPs = getSNPsFromDownstrNew(downstrFile)
        for markerAllele in downstrSNPs:
            w.write("\t".join([sampleName, markerAllele[0], markerAllele[1], ydnaSourceFile]) + "\n")
        
    w.close()

def putEscapeQuotes(string):
    return "\"" + string + "\""

def createUpdateScript(addAllelesFile, resultSummaryFileIn, resultSummaryFileOut, sampleName, mtDNASourceFile, yDNASourceFile, mtDNAFile, phyloEqFile, downstrFile, mtDNATag, phyloEqTag, downstrTag, novelSNPsTag, updateScriptFile):
    thisFileName = "getMTDNADifferences.py"
    with open(updateScriptFile, "w") as w:
        args = ["python3", thisFileName, "-update", resultSummaryFileIn, resultSummaryFileOut, mtDNAFile, phyloEqFile, downstrFile, putEscapeQuotes(mtDNATag), putEscapeQuotes(phyloEqTag), putEscapeQuotes(downstrTag), putEscapeQuotes(novelSNPsTag)]
        w.write(" ".join(args) + "\n")
        args = ["python3", thisFileName, "-addAlleles", addAllelesFile, mtDNAFile, sampleName, mtDNASourceFile, phyloEqFile, downstrFile, yDNASourceFile]
        w.write(" ".join(args) + "\n")
        args = ["cp", resultSummaryFileIn, resultSummaryFileIn + ".backup"]
        w.write(" ".join(args) + "\n")
        args = ["mv", resultSummaryFileOut, resultSummaryFileIn]
        w.write(" ".join(args) + "\n")
        args = ["python3", "create-mtDNA.py", sampleName]
        w.write(" ".join(args) + "\n")
    w.close()

#main part of the programm
#check for at least 4 parameter (len(sys.argv >= 5))
if len(sys.argv) > 4:

    #check for routine "process"
    if sys.argv[1] == "-process":
        vcfFile = sys.argv[2]
        haplogrepJar = sys.argv[3]
        outFile = sys.argv[4]
        if len(sys.argv) > 5:
            haplogrepOutputFile = sys.argv[5]
        process()
        
    if sys.argv[1] == "-update":
        resultsSummaryFileIn = sys.argv[2]
        resultsSummaryFileOut = sys.argv[3]
        mtDNAFile = sys.argv[4]
        phyloEqFile = sys.argv[5]
        downstrFile = sys.argv[6]
        mtdnaTagStr = sys.argv[7]
        phyloEqTagStr = sys.argv[8]
        downstrTagStr = sys.argv[9]
        novelSNPsTagStr = sys.argv[10]
        update(resultsSummaryFileIn, resultsSummaryFileOut, mtDNAFile, phyloEqFile, downstrFile, mtdnaTagStr, phyloEqTagStr, downstrTagStr, novelSNPsTagStr)
        
    if sys.argv[1] == "-addAlleles":
        addAllelesFile = sys.argv[2]
        mtdnaOutFile = sys.argv[3]
        sampleName = sys.argv[4]
        mtdnaSourceFile = sys.argv[5]
        phyloEqFile = sys.argv[6]
        downstrFile = sys.argv[7]
        ydnaSourceFile = sys.argv[8]
        addAlleles(mtdnaOutFile, addAllelesFile, sampleName, mtdnaSourceFile, phyloEqFile, downstrFile, ydnaSourceFile)
        
    if sys.argv[1] == "-createUpdateScript":
        addAllelesFile = sys.argv[2]
        resultSummaryFileIn = sys.argv[3]
        resultSummaryFileOut = sys.argv[4]
        sampleName = sys.argv[5]
        mtDNASourceFile = sys.argv[6]
        yDNASourceFile = sys.argv[7]
        mtDNAFile = sys.argv[8]
        phyloEqFile = sys.argv[9]
        downstrFile = sys.argv[10]
        mtDNATag = sys.argv[11]
        phyloEqTag = sys.argv[12]
        downstrTag = sys.argv[13]
        novelSNPsTag = sys.argv[14]
        updateScriptFile = sys.argv[15]
        createUpdateScript(addAllelesFile, resultSummaryFileIn, resultSummaryFileOut, sampleName, mtDNASourceFile, yDNASourceFile, mtDNAFile, phyloEqFile, downstrFile, mtDNATag, phyloEqTag, downstrTag, novelSNPsTag, updateScriptFile)
        
