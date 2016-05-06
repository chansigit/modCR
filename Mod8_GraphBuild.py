# coding=utf-8
from Aligners.Aligner import *
from Aligners.EndsFree import *
from Aligners.Pairwise import *
from Aligners.LocalAligner import *
from Aligners.Multiple import *
from Aligners.SimpleMultiple2 import *
from Utils.HashUtils import *
from Utils.FileUtils import *
from time import *
from random import *
from sys import *
import operator
import os

# **********************************************************************************************************************
# *****************************************   Variable Definition Section   ********************************************
# **********************************************************************************************************************
# Get user input parameters
if (len(argv) < 2 or len(argv) % 2 != 0):
    print "python CrisprDetector.py <INPUT-FILE-PATH> [options]"
    exit()

# Input file data
(DIRECTORY, FILE) = os.path.split(argv[1])
DIRECTORY += os.path.sep

dictionaryArguments = {argv[i]: argv[i+1] for i in range(2, len(argv), 2)}

# k-mer statistics parameters
KVAL      = int(dictionaryArguments["-k"]) if "-k" in dictionaryArguments else 23
threshold = int(dictionaryArguments["-t"]) if "-t" in dictionaryArguments else 45

# Overlap detection parameters
MIN_OVERLAP_FRACTION = float(dictionaryArguments["-o"]) if "-o" in dictionaryArguments else 0.3
SIDE_ZONE            = float(dictionaryArguments["-sd"]) if "-sd" in dictionaryArguments else 0.4

# CRISPR array parameters
MINIMUM_SPACERS       = int(dictionaryArguments["-n"]) if "-n" in dictionaryArguments else 2
MINIMUM_REPEAT_LENGTH = int(dictionaryArguments["-r"]) if "-r" in dictionaryArguments else 21
MAXIMUM_REPEAT_LENGTH = int(dictionaryArguments["-R"]) if "-R" in dictionaryArguments else 60
MINIMUM_SPACER_LENGTH = int(dictionaryArguments["-s"]) if "-s" in dictionaryArguments else 15
MAXIMUM_SPACER_LENGTH = int(dictionaryArguments["-S"]) if "-S" in dictionaryArguments else 75

# Sampling parameters
UPPER_BOUND_NUM_ISOLATED_NODES = int(dictionaryArguments["-bn"]) if "-bn" in dictionaryArguments else 50
UPPER_BOUND_NUM_EDGES_PER_NODE = int(dictionaryArguments["-be"]) if "-be" in dictionaryArguments else 250

# Overlap constants
SIDES = [PREFIX, SUFFIX]
MOD   = 2 ** 31 - 1

# CRISPR detection and clustering constants
MINIMUM_READS_CLUSTER_MERGE = 35 	# minimum number of shared containig reads between to mergeable k-kmers
SPACER_KMER = 5 					# length of signature of a spacer beginning/end
SPACER_SIM_FRAC = 0.75				# threshold for similarity of spacers (as a fraction of the length)
SPACER_SIM_WEIGHT = 3				# weight for base similarity for spacer alignment
CONSEC_DIFFS_FOR_SPACER = 2			# number of consecutive different bases searched for the determination of the spacer boundaries

# Consensus derivation constants
CONSENSUS_FRACTION_DOMINANT = 0.85  		# threshold for the fraction of reads having a dominant base in the consensus sequence
CONSENSUS_FRACTION_MULTIPLE_FIRST = 0.5		# threshold for the fraction of reads having the first non-dominant base in the consensus sequence
CONSENSUS_FRACTION_MULTIPLE_SECOND = 0.17	# threshold for the fraction of reads having the second non-dominant base in the consensus sequence
CONSEN_FRAC_MULT_MAX_THIRD = 0.05			# upper bound on the fraction of reads having a third non-dominant base in the consensus sequence

DEBUG = False

# Strand orientation constants
UNDEF = 0
REGULAR = 1
REVCOMP = 2
RC_DICT = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N', '{':'}', '}':'{', ' ':' '}

# **********************************************************************************************************************
# *****************************************   Function Definition Section   ********************************************
# **********************************************************************************************************************
# Get a list of overlapping reads and overlaps sizes for a given read
def GetOverlapDataForRead(readindex, read, data, minOverlap):
    overlaps = dict();
    hashValues = CalcHashValues(readindex, read, minOverlap, SUFFIX, None, MOD)
    for index in range(len(hashValues)):
        overlapSize = index + minOverlap
        hashValue = hashValues[index]
        if hashValue in data[overlapSize]:
            for overlapRead in data[overlapSize][hashValue][PREFIX]:
                overlaps[overlapRead] = overlapSize
    return overlaps

# Check if an edge is a spacer edge or not
def IsKmerInOverlap(left, right, kmer, overlapsSize, readLen):
    kval = len(kmer)
    leftIndices = []
    for index in range(len(left) - kval + 1):
        currkmer = left[index:index + kval]
        if (currkmer == kmer):
            leftIndices.append(index)
    rightIndices = []
    for index in range(len(right) - kval + 1):
        currkmer = right[index:index + kval]
        if (currkmer == kmer):
            rightIndices.append(index)
    for leftIndex in leftIndices:
        for rightIndex in rightIndices:
            if (leftIndex == rightIndex + (readLen - overlapsSize)):
                return True
    return False

# A binary search method
def binary_search(a, x, lo=0, hi=None):
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        midval = a[mid]
        if midval < x:
            lo = mid+1
        elif midval > x:
            hi = mid
        else:
            return mid
    return -1

# Get the offset of the spacers in the edge sequences
def GetSpacersOffsets(spacers):
    offsetRight = 0
    offsetLeft = 0
    minLen = min([len(spacer) for spacer in spacers])
    consecutiveDiffs = 0
    # Scan the edge sequences from the beginning
    for index in range(minLen):
        isDiff = False
        base = spacers[0][index]
        for sIndex in range(1, len(spacers)):
            if (spacers[sIndex][index] != base):
                consecutiveDiffs += 1
                isDiff = True
                break
        if (not isDiff):
            consecutiveDiffs = 0
        if (consecutiveDiffs == CONSEC_DIFFS_FOR_SPACER):
            offsetRight = index + 1
            break
    consecutiveDiffs = 0
    # Scan the edge sequences from the end
    for index in range(1, minLen):
        isDiff = False
        base = spacers[0][len(spacers[0]) - index]
        for sIndex in range(1, len(spacers)):
            if (spacers[sIndex][len(spacers[sIndex]) - index] != base):
                consecutiveDiffs += 1
                isDiff = True
                break
        if (not isDiff):
            consecutiveDiffs = 0
        if (consecutiveDiffs == CONSEC_DIFFS_FOR_SPACER):
            offsetLeft = index
            break
    return (offsetRight, offsetLeft)

# Generic method for spacer sdjacency check
def CompareSpacersLinkage(first, second, verificationIndexStartSpacer, verificationIndexEndSpacer, offsetLeft, offsetRight, spacerKval):
    start = second[(offsetRight - 1) : (offsetRight - 1) + spacerKval]
    end = first[-offsetLeft + 1 - spacerKval : -offsetLeft + 1]
    if (start in verificationIndexStartSpacer and end in verificationIndexEndSpacer and verificationIndexStartSpacer[start] == verificationIndexEndSpacer[end]):
        return True
    return False

# Cache all reads by the offsets of the spacers
# Caching is done by short kmers from the beginning and end of the spacers
def GetSpacerKmerCache(offsetRight, offsetLeft, kmer, spacerKval, readIndices, reads, verificatioIndexStartSpacer, verificatioIndexEndSpacer):
    readLength = len(reads[0])
    for readIndex in readIndices:
        read = reads[readIndex]
        kIndex = read.find(kmer)
        rightShift = KVAL + (offsetRight - 1) + spacerKval
        leftShift = spacerKval + offsetLeft - 1
        if (kIndex >= leftShift and kIndex <= readLength - rightShift):
            verificatioIndexStartSpacer[read[kIndex + rightShift - spacerKval : kIndex + rightShift]] = readIndex
            verificatioIndexEndSpacer[read[kIndex - leftShift : kIndex - leftShift + spacerKval]] = readIndex

# Conduct an inital check of adjacency between spacers (all pairs)
def CompareInitialSpacers(spacers, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval):
    matches = 0
    for firstIndex in range(len(spacers)):
        firstSpacer = spacers[firstIndex]
        for secondIndex in range(firstIndex + 1, len(spacers)):
            secondSpacer = spacers[secondIndex]
            if (CompareSpacersLinkage(firstSpacer, secondSpacer, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval)):
                matches += 1
            elif (CompareSpacersLinkage(secondSpacer, firstSpacer, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval)):
                matches += 1
            if (matches == max(MINIMUM_SPACERS - 1,1)):
                return matches
    return matches

# Conduct an adjacency check of a new spacer with previously detected ones
def CompareNewSpacer(nspacer, spacers, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval, matchesLeft):
    matches = 0
    for index in range(1,len(spacers)):
        otherspacer = spacers[index]
        if (CompareSpacersLinkage(nspacer, otherspacer, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval)):
            matches += 1
        elif (CompareSpacersLinkage(otherspacer, nspacer, verificatioIndexStartSpacer, verificatioIndexEndSpacer, offsetLeft, offsetRight, spacerKval)):
            matches += 1
        if (matches == matchesLeft):
            return matches
    return matches

# Analyze a kmer and decide whether it belongs to a CRISPR repeat
def AnalyzeKmer(kmer, readIndices, reads, overlapsData, cacheOverlapsData, minOverlap):
    readLen = len(reads[0])
    readNumber = len(readIndices)
    startReads = []
    endReads = []
    # Go over all reads and mark reads that contain the kmer towards the beginning or end
    for readIndex in readIndices:
        read = reads[readIndex]
        position = read.find(kmer)
        if (position > -1 and position < readLen * SIDE_ZONE):
            startReads.append(readIndex)
        if (position > readLen * (1 - SIDE_ZONE)):
            endReads.append(readIndex)
    # Variables definitions
    spacers = []
    duplicateSpacers = []
    (offsetRight, offsetLeft) = (None, None)
    verificationIndexStartSpacer = dict()
    verificationIndexEndSpacer = dict()
    currentStarts = []
    currentEnds = []
    badStarters = 0
    goodStarters = 0
    endReads.sort()

    # Go over reads that contain the kmer at the beginning and find suitable overlaps
    for startIndex in startReads:
        # Check if upper bound on number of isolated reads has been reached
        if (badStarters > UPPER_BOUND_NUM_ISOLATED_NODES and goodStarters < UPPER_BOUND_NUM_ISOLATED_NODES):
            return []
        # Define a flag indicating whether the current reads is isolated or not
        isBadStarter = True
        # Obtain overlap data for current read
        if (startIndex not in cacheOverlapsData):
            cacheOverlapsData[startIndex] = GetOverlapDataForRead(startIndex, reads[startIndex], overlapsData, minOverlap)
        overlaps = cacheOverlapsData[startIndex]
        # Go over all reads overlapping current read to the right
        matches = 0
        left = reads[startIndex]
        kmerIndexInStart = left.find(kmer)
        assert (kmerIndexInStart >  -1)
        enders = 0
        for endIndex in overlaps:
            overlapSize = overlaps[endIndex]
            # If the overlap is too large then the overlap is not related to a spacer edge
            if (overlapSize >= readLen - kmerIndexInStart):
                continue
            # Check if upper bound on the number of edges per node has been reached
            if (enders > UPPER_BOUND_NUM_EDGES_PER_NODE):
                break
            # Update the edges counter
            enders += 1
            # Check that the overlapping read contains the kmer towards the end
            if (binary_search(endReads, endIndex) > -1):
                if (startIndex == endIndex):
                    continue
                right = reads[endIndex]
                # Verify that the overlap corresponds to a spacer edge
                kmerInOverlap = IsKmerInOverlap(left, right, kmer, overlaps[endIndex], readLen)
                if (not kmerInOverlap):
                    # If there is at least one overlap corresponding to a spacer edge - this is a good starter
                    # Extract edge sequence
                    sequence = left + right[overlapSize:]
                    spacer = sequence[left.find(kmer) + len(kmer) : readLen - overlapSize + right.find(kmer)]
                    if (len(spacer) < MINIMUM_SPACER_LENGTH or len(spacer) > MAXIMUM_SPACER_LENGTH):
                        continue
                    isBadStarter = False
                    # Check that the edge sequence is not already known (fast check)
                    if spacer not in spacers and spacer not in duplicateSpacers and kmer not in spacer:
                        shouldAdd = True
                        # Verify that the edge sequence is new (pairwise alignment)
                        for knownSpacer in spacers:
                            aligner = Pairwise(spacer, knownSpacer)
                            alignment = aligner.align()
                            if (alignment.score > len(spacer) * SPACER_SIM_FRAC * SPACER_SIM_WEIGHT - len(spacer) * (1 - SPACER_SIM_FRAC)):
                                shouldAdd = False
                                duplicateSpacers.append(spacer)
                                break
                        # Add spacer if indeed not already known
                        if (shouldAdd):
                            spacers.append(spacer)
                            # Get spacer offsets (if enough spacers were already found)
                            if (len(spacers) == MINIMUM_SPACERS):
                                (offsetRight, offsetLeft) = GetSpacersOffsets(spacers)
                                if (offsetRight == None or offsetLeft == None):
                                    print ("Error in identifying spacer k-mer!!!")
                                    assert(False)
                                # Cache reads by spacer offsets
                                GetSpacerKmerCache(offsetRight, offsetLeft, kmer, SPACER_KMER, readIndices, reads, verificationIndexStartSpacer, verificationIndexEndSpacer)
                               # Find initial number of consecutive spacers
                                matches = CompareInitialSpacers(spacers, verificationIndexStartSpacer, verificationIndexEndSpacer, offsetLeft, offsetRight, SPACER_KMER)
                            # If enough consecutive spacers were found then the kmer belongs to a CRISPR repeat
                            if (matches >= max(MINIMUM_SPACERS - 1,1)):
                                return spacers
                            # Continue searching for adjacency between spacers (if new ones are found)
                            if (len(spacers) > MINIMUM_SPACERS):
                                matchesLeft = max(MINIMUM_SPACERS - 1,1) - matches
                                matches += CompareNewSpacer(spacer, spacers, verificationIndexStartSpacer, verificationIndexEndSpacer, offsetLeft, offsetRight, SPACER_KMER, matchesLeft)
                                if (matches >= max(MINIMUM_SPACERS - 1,1)):
                                    return spacers
        # Update counters of isolated and non-isolated nodes
        if (isBadStarter):
    	    badStarters += 1
        else:
            goodStarters += 1
    spacers = []
    return spacers


# **********************************************************************************************************************
# *********************************************   Computation Section   ************************************************
# **********************************************************************************************************************
# python Mod8_GraphBuild.py c:\data\chromosome.fasta -k 20 -t 40 -retdir c:\data -log c:\data\ -tname may1test
# -rfst C:\data\Mod4_RefinedStats_may05test_bcc77d80-11fc-11e6-bac8-ec55f98094e4.json
# -rfrd C:\data\Mod5_RefinedReads_may05test_e91ab04f-11fc-11e6-870d-ec55f98094e4.json
# -hash C:\data\Mod6_HashData_may05test_bd9fa840-12b0-11e6-b637-ec55f98094e4.json
# -pair C:\data\Mod7_MatchedPair_may05test_55e83c21-11fd-11e6-bbf5-ec55f98094e4.json
# -ovlp C:\data\Mod3_minOverlap_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json
# -rlrd C:\data\Mod3_relevantReadsNumber_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json
# -rlen C:\data\Mod3_readLen_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json
# python Mod8_GraphBuild.py c:\data\chromosome.fasta -k 20 -t 40 -retdir c:\data -log c:\data\ -tname may1test -rfst C:\data\Mod4_RefinedStats_may05test_bcc77d80-11fc-11e6-bac8-ec55f98094e4.json  -rfrd C:\data\Mod5_RefinedReads_may05test_e91ab04f-11fc-11e6-870d-ec55f98094e4.json -hash C:\data\Mod6_HashData_may05test_bd9fa840-12b0-11e6-b637-ec55f98094e4.json -pair C:\data\Mod7_MatchedPair_may05test_55e83c21-11fd-11e6-bbf5-ec55f98094e4.json -ovlp C:\data\Mod3_minOverlap_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json -rlrd C:\data\Mod3_relevantReadsNumber_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json -rlen C:\data\Mod3_readLen_may05test_594629a1-11fc-11e6-8b62-ec55f98094e4.json

import json
import uuid
import logging
# 定义辅助参数
RESULTDIR= str(dictionaryArguments["-retdir"]) if "-retdir" in dictionaryArguments else "~/"
LOGDIR   = str(dictionaryArguments["-log"])    if "-log"    in dictionaryArguments else RESULTDIR
TASKNAME = str(dictionaryArguments["-tname"])  if "-tname"  in dictionaryArguments else "default"
refinedStatsFile=str(dictionaryArguments["-rfst"]) if "-rfst"  in dictionaryArguments else "refinedstat.json"
refinedReadsFile=str(dictionaryArguments["-rfrd"]) if "-rfrd"  in dictionaryArguments else "refinedread.json"
pairsFile= str(dictionaryArguments["-pair"])   if "-pair"  in dictionaryArguments else "pairs.json"
hashFile = str(dictionaryArguments["-hash"])   if "-hash"  in dictionaryArguments else "hash.json"
minOverlapFile = str(dictionaryArguments["-ovlp"])   if "-ovlp"  in dictionaryArguments else "minoverlap.json"
relevantReadsNumberFile = str(dictionaryArguments["-rlrd"])   if "-rlrd"  in dictionaryArguments else "relevantRead.json"
readLenFile = str(dictionaryArguments["-rlen"])   if "-rlen"  in dictionaryArguments else "readLen.json"
# 生成日志文件名和输出文件名
uuidstr=str(uuid.uuid1())
Ret_Filename = "Mod8_GraphBuild_"+TASKNAME+"_"+uuidstr+".json"  #
Ret_Filename = os.path.join(RESULTDIR, Ret_Filename)            #
Log_Filename = "Mod8_GraphBuild_"+TASKNAME+"_"+uuidstr+".log"
Log_Filename = os.path.join(LOGDIR, Log_Filename)

# 写入任务信息
logging.basicConfig(level = logging.DEBUG, datefmt = '%a, %d %b %Y %H:%M:%S', filename = Log_Filename, filemode = 'w',
                    format = '%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
logtitle = "Module GraphBuild - "+TASKNAME
logging.info(logtitle)
logging.info("Parameter K=" + str(KVAL) + ", T=" +str(threshold))
logging.info("Logging in "+str(Log_Filename))
logging.info("Results in "+str(Ret_Filename))

# 导入数据
logging.info("File Loading begins")
t1_load = time()
refinedStats = json.load(open(refinedStatsFile, 'r'))
refinedReads = json.load(open(refinedReadsFile, 'r'))
pairs        = json.load(open(pairsFile, 'r'))
data         = json.load(open(hashFile,  'r'))
minOverlap          = json.load(open(minOverlapFile,'r'))
relevantReadsNumber = json.load(open(relevantReadsNumberFile, 'r'))
readLen             = json.load(open(readLenFile, 'r'))
t2_load = time()
logging.info("File Loading Finished, taking " + str(t2_load-t1_load) + " seconds")

# 计算
logging.info("Graph Building begins")
t1_GraphBuild = time()
# -------------------------------------------
overlapsCache = dict()
crisprKmers   = []
kindex        = 0
gapt          = 0
bapt          = 0
# -------------------------------------------
startTimeKmerAnalysis = time()
for (kmer, rckmer) in pairs.items():
    readIndices = refinedStats[kmer]
    rcReadIndices = []
    if rckmer != None:
        rcReadIndices = refinedStats[rckmer]

    allIndices = []
    for regularIndex in readIndices:
        allIndices.append(regularIndex)
    for rcIndex in rcReadIndices:
        allIndices.append(relevantReadsNumber + rcIndex)

    apt1 = time()
    print "kmer=", kmer
    spacers = AnalyzeKmer(kmer, allIndices, refinedReads, data, overlapsCache, minOverlap)
    print spacers
    apt2 = time()
    kindex += 1

    if (len(spacers) > 1):
        crisprKmers.append((kmer, rckmer))
        gapt += apt2 - apt1
    else:
        bapt += apt2 - apt1
endTimeKmerAnalysis = time()
logging.info("Analyzing kmer takes ",endTimeKmerAnalysis - startTimeKmerAnalysis, "seconds")
# -------------------------------------------

# -------------------------------------------
t2_GraphBuild = time()
logging.info("Graph Building ends, taking "+ str(t2_GraphBuild-t1_GraphBuild) +" seconds")



# 持久化结果
logging.info("Generating JSON file begins")
fp = open(Ret_Filename,"w")
t1_json =time()

#json.dump(xxx, fp)

t2_json =time()
fp.close()
logging.info("Generating JSON file ends, taking "+ str(t2_json-t1_json) +" seconds")
