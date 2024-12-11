# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:40:12 2020

@author: hunte
"""

import json
import sys
import vcf

#yfullFile = "C:\YSEQ\FASTQ to hg38 Mapping\scripts\latest_YFull_YTree.json"
clade = "J-FGC58561"

if len(sys.argv) > 5:
    yfullFile = sys.argv[1]
    clade = sys.argv[2].replace("*","")
    vcfFile = sys.argv[3]
    positivesFile = sys.argv[4]
    negativesFile = sys.argv[5]
    phyloEquivFile = sys.argv[6]
    downstreamFile = sys.argv[7]
    idxStatsTSVFile = sys.argv[8]

    
obj=json.load(open(yfullFile))

phyloequivalentSNPs = []

def findClade(clade, node):
    if "id" not in node or node["id"] != clade:
        if "children" in node:
            for child in node["children"]:
                res = findClade(clade, child)
                if res != None:
                    return res
        return None
    return node

def snpStringToArray(s):
    snps = s.split(", ")
    if "" in snps:
        snps.remove("")
    return snps

def getDownstreamSNPs(node):
    snps = {}
    if "children" in node:
        for child in node["children"]:
            if "*" not in child["id"]:
                snps[child["id"]] = set(snpStringToArray(child["snps"]))
    return snps

foundClade = findClade(clade, obj)
phyloequivalentSNPs = snpStringToArray(foundClade['snps'])
downstreamBranches = getDownstreamSNPs(foundClade)
downstream = set([])
for branch in downstreamBranches:
    downstream = downstream.union(downstreamBranches[branch])

if len(downstream) == 0:
    downstream = set(["None"])

downstreamDict = {}
for down in downstream:
    downstreamDict[down] = []
    for spl in down.split("/"):
        downstreamDict[down].append(spl)

phyloEquivDict = {}
for phylo in phyloequivalentSNPs:
    phyloEquivDict[phylo] = []
    for spl in phylo.split("/"):
        phyloEquivDict[phylo].append(spl)
    
#print("\t".join(["PHYLOEQUIVALENT_SNPS","DOWNSTREAM_SNPS"]))
#print("\t".join([", ".join(phyloequivalentSNPs), ", ".join(downstream)]))

def getEquivalent(inputPossiblePhyloEqSNPs, thedict):
    for snp in inputPossiblePhyloEqSNPs:
        for dictkey in thedict:
            if snp in thedict[dictkey]:
                return dictkey
    return None


phyloEquivFound = {}
for phy in phyloEquivDict:
    phyloEquivFound[phy] = "not found"

downEquivFound = {}
for down in downstreamDict:
    downEquivFound[down] = "not found"

positives = set([])
negatives = set([])


def isMale():

    chrXcov = 0
    chrYcov = 1
    with open(idxStatsTSVFile, "r") as r:
        lines = r.readlines()
       
        for line in lines:
            split = line.split("\t")
            thechr = split[0]
            if thechr == "chrY":
                chrYcov = float(split[2]) / float(split[1])
            if thechr == "chrX":
                chrXcov = float(split[2]) / float(split[1])
    r.close()
    if chrXcov / chrYcov > 5:
        return False
    return True

sampleName = idxStatsTSVFile.split(".idxstats")[0]
def getAllele(record):
    basesString = record.genotype(sampleName).gt_bases
    if basesString:
        basesSplits = basesString.split("/")
        if len(basesSplits) > 1:
            call1 = basesSplits[0]
            call2 = basesSplits[1]
            if call1 == call2:
                return call1
            return "mixed"
        else:
            return "not parseable"
    return "not parseable"
        
        
if isMale():    
    vcf_reader = vcf.Reader(filename=vcfFile)
    vcf_reader.fetch('chrY')

    for record in vcf_reader:
        theids = record.ID.split(",")
        phyloEquiv = getEquivalent(theids, phyloEquivDict)
        if phyloEquiv:
            phyloEquivFound[phyloEquiv] = getAllele(record)
        down = getEquivalent(theids, downstreamDict)
        if down:
            downEquivFound[down] = getAllele(record)
    
    def readSNPsFromFile(file):
        with open(file, "r") as r:
            return r.readline().split(",")
        
    positives = set(readSNPsFromFile(positivesFile))
    negatives = set(readSNPsFromFile(negativesFile))
    if '' in negatives:
        negatives.remove('')
    if '' in positives:
        positives.remove('')
    allRelevant = set([])
    for phy in phyloEquivDict:
        for pe in phyloEquivDict[phy]:
            allRelevant.add(pe)
    for d in downstreamDict:
        for de in downstreamDict[d]:
            allRelevant.add(de)

    positives = positives.intersection(allRelevant)
    negatives = negatives.intersection(allRelevant)
else:
    phyloEquivFound = {}
    downstreamBranches = {}

with open(phyloEquivFile, "w") as w:
    w.write("\t".join(["position","mutation","call"]) + "\n")
    for phy in phyloEquivFound:
        call = "no call"
        if any(p in positives for p in phyloEquivDict[phy]):
            call = "+"
        if any(p in negatives for p in phyloEquivDict[phy]):
            call = "-"
        mutation = "unknown"
        if phyloEquivFound[phy] not in ["not found", "mixed", "not parseable"]:
            mutation = phyloEquivFound[phy]
        w.write("\t".join([phy, mutation, call]) + "\n")
w.close()


with open(downstreamFile, "w") as w:
    w.write("\t".join(["branch","position","mutation","call"]) + "\n")
    for branch in downstreamBranches:
        for down in downstreamBranches[branch]:
            call = "no call"
            if any(d in positives for d in downstreamDict[down]):
                call = "+"
            if any(d in negatives for d in downstreamDict[down]):
                call = "-"
            mutation = "unknown"
            if downEquivFound[down] not in ["not found", "mixed", "not parseable"]:
                mutation = downEquivFound[down]
            w.write("\t".join([branch, down, mutation, call]) + "\n")
w.close()