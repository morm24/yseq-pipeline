# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:23:13 2018

@author: hunte
"""

"""
python cladeFinder.py YFull_YTree_v6.02__20180402.tree.json pos negatives output

create negatives
bcftools filter -O z -i '(GT=="1/1" && AA==ALT) || (GT=="0/0" && AA=REF)' chrY_cleaned_1_hg38.vcf | bcftools query -f '%ID,' > negatives

create positives
bcftools query -f '%ID,' chrY_derived_1_hg38.vcf.gz > pos

"""


import json
import sys



toIgnore = ["PF129", "Z2533", "S6868", "BY2285", "YP2229", "YP2250", "YP2228", "YP2129", "YP1838", "YP1807", "YP1841", "YP1740", "Y17293", "PF6234",  "L132.2"]
treeFile = "/genomes/0/tree/trees/yfull/latest_YFull_YTree.json"
positivesFile = "C://Users//hunte//Documents//YSeq//pos"
negativesFile = "C://Users//hunte//Documents//YSeq//negatives"
outputFile = "C://Users//hunte//Documents//YSeq//cladeDetermination.csv"


if len(sys.argv) > 1:
    treeFile = sys.argv[1]
    positivesFile = sys.argv[2]
    negativesFile = sys.argv[3]
    outputFile = sys.argv[4]
 
def parseFile(file):
    fr = open(file, 'r')
    return fr.readline().replace("\n","").split(",")

positives = set(parseFile(positivesFile))
negatives = set(parseFile(negativesFile))
if "" in positives:
    positives.remove("")
if "" in negatives:
    negatives.remove("")
for ign in toIgnore:
    if ign in positives:
        positives.remove(ign)
    if ign in negatives:
        negatives.remove(ign)

hierarchy = {}
snps = {}

snpToBranch = {}

def parseTreeJSON(fil):
    root = json.load(open(fil))
    recurseTreeJson(root, hierarchy, snps)
    return (root["id"], hierarchy, snps)

def parseSNPsString(branch, snpsString):
    thesnps = []
    for snps in snpsString.split(", "):
        for snp in snps.split("/"):
            thesnps.append(snp)
            if snp not in snpToBranch:
                snpToBranch[snp] = [branch]
            else:
                snpToBranch[snp].append(branch)
    return thesnps

def getPossibleTerminals():
    possible = set([])
    for pos in positives:
        if pos in snpToBranch:
            for branch in snpToBranch[pos]:
                possible.add(branch)
    return possible
            
def recurseTreeJson(node, hierarchy, snps):
    if "children" in node:
        for child in node["children"]:
            hierarchy[child["id"]] = node["id"]
            snps[child["id"]] = parseSNPsString(child["id"],child["snps"])
            recurseTreeJson(child, hierarchy, snps)
                

def getChildren(clade, childParents):
    children = []
    for child in childParents:
        if childParents[child] == clade:
            children.append(child)

    return children

#positives = ["M241","Z1043", "M269", "Z1297", "PH1553", "Z8424", "Y27522"]
#negatives = ["L283", "M231"]

def isInChildrenThisLevel(clade, positives, childParents):
    children = getChildren(clade, childParents)
    inChildren = []
    for child in children:
        if any(snp in positives for snp in snps[child]):
            inChildren.append(child)
    return inChildren

def recurseDownTreeUntilFirstHits(clade, positives, childParents):
    posChildrenThisLevel = isInChildrenThisLevel(clade, positives, childParents)
    for child in getChildren(clade, childParents):
        if child not in posChildrenThisLevel:
            childResult = recurseDownTreeUntilFirstHits(child, positives, childParents)
            for cres in childResult:
                posChildrenThisLevel.append(cres)
    return posChildrenThisLevel

solutions = []

def removeDuplicates(arr): 

    n = len(arr)
    # Return, if array is  
    # empty or contains 
    # a single element 
    if n == 0 or n == 1: 
        return n 
  
    temp = list(range(n)) 
  
    # Start traversing elements 
    j = 0; 
    for i in range(0, n-1): 
  
        # If current element is 
        # not equal to next 
        # element then store that 
        # current element 
        if arr[i] != arr[i+1]: 
            temp[j] = arr[i] 
            j += 1
  
    # Store the last element 
    # as whether it is unique 
    # or repeated, it hasn't 
    # stored previously 
    temp[j] = arr[n-1] 
    j += 1
      
    # Modify original array 
    for i in range(0, j): 
        arr[i] = temp[i] 
  
    return arr

def refineHitsRecursively(sequences, positives, childParents):
    for sequence in sequences:
        refinedResults = recurseDownTreeUntilFirstHits(sequence[-1], positives, childParents)
        if len(refinedResults) == 0:
            solutions.append(sequence)
        else:
            print(sequence, refinedResults)
            for refRes in refinedResults:
                #print(sequence, refRes)
                seqCopy = sequence[:]
                seqCopy.append(refRes)
                refineHitsRecursively([seqCopy], positives, childParents)               

def recurseDownTree(positives, childParents):
    sequences = recurseDownTreeUntilFirstHits("", positives, childParents)
    newSequences = []
    for sequence in sequences:
        newSequences.append([sequence])
    refineHitsRecursively(newSequences, positives, childParents)

def getTotalSequence(clade, hierarchy):
    sequence = [clade]
    thisClade = clade
    while thisClade in hierarchy:
        thisClade = hierarchy[thisClade]
        sequence.append(thisClade)
    return sequence[:-1]
    
def getScore(sequence, totalSequence):
    return float(len(sequence)) / len(totalSequence)

def printSolutions(solutions):
    for solution in solutions:
        print(" ".join(solution), getScore(solution))

def getConflicts(sequence, negatives, hierarchy):
    conflictingNegatives = []
    hgs = []
    for hg in sequence:
        if any(snp in negatives for snp in snps[hg]):
            conflictingNegativeSnps = ""
            for snp in snps[hg]:
                if snp in negatives:
                    conflictingNegativeSnps += " " + snp
            conflictingNegatives.append(hg + " @" + conflictingNegativeSnps + ";")
            hgs.append(hg)
    return conflictingNegatives, hgs

def getWarnings(sequence, negatives, hierarchy):
    messages = []
    (conflicts, _) = getConflicts(sequence, negatives, hierarchy)
    for conflict in conflicts:
        messages.append(" " + conflict)
    return messages

def getWarningsConf(conflicts):
    messages =[]
    for conflict in conflicts:
        messages.append(" " + conflict)
    return messages

def unconflictedPercent(sequence, conflicts):
    return (float(len(sequence)) - len(conflicts)) / len(sequence)

import numpy as np
def getPathScores(fullSequence, confirmed, negatives, positives, conflicts):
    scores = []
    last = fullSequence[-1]
    weights = 0
    for thing in fullSequence:
        weight = len(snps[thing])
        if thing in confirmed:
            negs = len(negatives.intersection(set(snps[thing])))
            poses = len(positives.intersection(set(snps[thing])))
            if last == thing:
                scores.append(1.0 * weight)
                weights += weight
            else:
                scores.append(weight * (-1.0 + 2.0 * float(poses) / float(poses + negs)))
                weights += weight
        else:
            if thing in conflicts:
                negs = len(negatives.intersection(set(snps[thing])))
                scores.append(weight * -1 * negs)
                weights += weight
    return np.divide(scores, weights)

def isBasal(clade, negatives, positives, hierarchy):
    basal = False
    children = getChildren(clade, hierarchy)
    if len(children) > 0:
        basal = True
        for child in children:
            isNeg = len(negatives.intersection(set(snps[child]))) > 0
            isPos = len(positives.intersection(set(snps[child]))) > 0
            if isNeg and not isPos:
                basal = basal and True
    return basal        

def writeFile(scoredSolutions):
    w = open(outputFile,"w+")
    w.write("#Haplogroup\tPath\tPurity\tDepth\tScore\tConflicting Negative SNPs\n")
    for scoredSolution in scoredSolutions:
        w.write(scoredSolution[1] + "\t" + " > ".join(scoredSolution[0]) + "\t" + str(scoredSolution[2]) + "\t" + str(scoredSolution[3]) + "\t" + str(scoredSolution[4]) + "\t" + " ".join(scoredSolution[5]))
        w.write("\n")
    w.close()

def getPositive(scores):
    totscore = 0
    for score in scores:
        if score > 0:
            totscore += 1
    return totscore

from operator import itemgetter
def getRankedSolutions(positives, negatives, hierarchy):
    scoredSolutions = []
    uniqueSolutions = getPossibleTerminals()
    for solution in uniqueSolutions:
        totalSequence = getTotalSequence(solution, hierarchy)
        totalSequence.reverse()
        (conflicts, conflictingHGs) = getConflicts(totalSequence, negatives, hierarchy)
        scores = getPathScores(totalSequence, uniqueSolutions, negatives, positives, conflictingHGs)
        clade = solution
        if isBasal(clade, negatives, positives, hierarchy):
            clade = clade + "*"
        scoredSolutions.append([totalSequence, clade, np.sum(scores), len(scores), getPositive(scores) * np.sum(scores), getWarningsConf(conflicts)])
        print(totalSequence, scores, str(np.average(scores)))
            
    scoredSolutions = sorted(scoredSolutions, key=itemgetter(4), reverse=True)
    
    scoredSolutions[0]
    return scoredSolutions
        
a = parseTreeJSON(treeFile)
hierarchy = a[1]
snps = a[2]
print(hierarchy["A00c"])
print(snps["A00"])
print(getChildren("",hierarchy))
print(snpToBranch["PH1080"])
print(getPossibleTerminals())
b = getRankedSolutions(positives, negatives, hierarchy)
writeFile(b)
