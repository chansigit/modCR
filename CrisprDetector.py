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

# Get user input parameters
if (len(argv) < 2 or len(argv) % 2 != 0):
    print "python CrisprDetector.py <INPUT-FILE-PATH> [options]"
    exit()

# Input file data	
(DIRECTORY, FILE) = os.path.split(argv[1])
DIRECTORY += os.path.sep

dictionaryArguments = {argv[i]: argv[i+1] for i in range(2, len(argv), 2)}

# k-mer statistics parameters
KVAL = int(dictionaryArguments["-k"]) if "-k" in dictionaryArguments else 23
threshold = int(dictionaryArguments["-t"]) if "-t" in dictionaryArguments else 45

# Overlap detection parameters
MIN_OVERLAP_FRACTION = float(dictionaryArguments["-o"]) if "-o" in dictionaryArguments else 0.3
SIDE_ZONE = float(dictionaryArguments["-sd"]) if "-sd" in dictionaryArguments else 0.4

# CRISPR array parameters
MINIMUM_SPACERS = int(dictionaryArguments["-n"]) if "-n" in dictionaryArguments else 2
MINIMUM_REPEAT_LENGTH = int(dictionaryArguments["-r"]) if "-r" in dictionaryArguments else 21
MAXIMUM_REPEAT_LENGTH = int(dictionaryArguments["-R"]) if "-R" in dictionaryArguments else 60
MINIMUM_SPACER_LENGTH = int(dictionaryArguments["-s"]) if "-s" in dictionaryArguments else 15
MAXIMUM_SPACER_LENGTH = int(dictionaryArguments["-S"]) if "-S" in dictionaryArguments else 75

# Sampling parameters
UPPER_BOUND_NUM_ISOLATED_NODES = int(dictionaryArguments["-bn"]) if "-bn" in dictionaryArguments else 50
UPPER_BOUND_NUM_EDGES_PER_NODE = int(dictionaryArguments["-be"]) if "-be" in dictionaryArguments else 250

# Overlap constants
SIDES = [PREFIX, SUFFIX]
MOD = 2 ** 31 - 1

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

# Get the length of the reads (assuming all have the same length)
def GetReadsLength(directory, fileName):
    handler = open(directory + fileName, "r")
    line = handler.readline().strip()
    while line.startswith(">"):
        line = handler.readline().strip()
    return len(line);
	

# Obtain frequency statistics on kmers
# kmers that appeat more than once in a read are currently discarded
def GetKmerStats(reads, kval):
    stats = dict()
    readindex = 0
    for read in reads:
        kmersInRead = set()
        if (readindex % 100000 == 0):
            stdout.write("\r Getting k-mer statistics. processed: %d" % readindex)
            stdout.flush()
        for index in range(len(read) - kval + 1):
            kmer = read[index:index + kval]
            # Check that kmer appears once in the read
            if (kmer in kmersInRead):
                if (kmer in stats):
                    del stats[kmer]
            else:
                kmersInRead.add(kmer)
                stats[kmer] = stats.get(kmer, 0) + 1
        readindex += 1
        del kmersInRead
    stdout.write("\r Getting k-mer statistics. processed: %d" % readindex)
    stdout.flush()
    stdout.write("\n")
    return stats


# Reverse complement a sequence
def RC(seq):
    return ''.join([RC_DICT[base] for base in seq[::-1]])

	
# Refine kmer statistics object to hold data only on frequent kmers
# In addition, the value is changed to hold list for read indices
def GetRefineStats(stats, threshold):
    refinedStats = dict()
    for kmer in stats:
        amount = stats[kmer]
        rc = RC(kmer)
        if rc in stats:
            amount += stats[rc]
        if (amount >= threshold):
            refinedStats[kmer] = list()
            refinedStats[rc] = list()
    return refinedStats


# Get all reads that contain some frequent kmer
def GetRelevantReads(directory, fileName, kval, data):
    signatures = set()
    refinedreads = []
    readindex = 0
    reads = readFileReadsThin(directory + fileName)
    for read in reads:
        if (readindex % 100000 == 0):
            stdout.write("\r Getting relevant reads. processed: %d" % readindex)
            stdout.flush()
        hashValue = hash(read)
        if hashValue in signatures:
            continue
        signatures.add(hashValue)		
        for index in range(len(read) - kval + 1):
            kmer = read[index:index + kval]
            if kmer in data:
                refinedreads.append(read)
                break
        readindex += 1
    stdout.write("\r Getting relevant reads. processed: %d" % readindex)
    stdout.flush()
    stdout.write("\n")
    return refinedreads;


# Obtain kmer index data for overlap queries
def ComputeHashValues(reads, minOverlap, mod):
    data = dict()
    readindex = 0
    for read in reads:
        if (readindex % 100000 == 0):
            stdout.write("\r Computing hash values. processed: %d" % readindex)
            stdout.flush()
        for side in [PREFIX]:
            hashValues = CalcHashValues(readindex, read, minOverlap, side, None, MOD)
            for index in range(len(hashValues)):
                sampleIndex = index + minOverlap
                hashValue = hashValues[index]
                data[sampleIndex] = data.get(sampleIndex, dict())
                data[sampleIndex][hashValue] = data[sampleIndex].get(hashValue, ([],[]))
                data[sampleIndex][hashValue][side].append(readindex)            
        readindex += 1
    stdout.write("\r Getting relevant reads. processed: %d" % readindex)
    stdout.flush()
    stdout.write("\n")
    return data
   
 
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

	
# Duplicate a set of reads by appending the reverse complement reads
def DuplicateReads(reads):
    originals = len(reads)
    for readIndex in range(len(reads)):
        reads.append(RC(reads[readIndex]))
    return reads


# Find pairs of reads and their reverse complements
def MatchKmers(stats):
    pairs = dict()
    for kmer in stats:
        rc = RC(kmer)
        if (rc in stats):
            if rc not in pairs:
                pairs[kmer] = rc
        else:
            pairs[kmer] = None
    return pairs


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
                                    assert(false)
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


# Generic method for spacer sdjacency check
def CompareSpacersLinkage(first, second, verificationIndexStartSpacer, verificationIndexEndSpacer, offsetLeft, offsetRight, spacerKval):
    start = second[(offsetRight - 1) : (offsetRight - 1) + spacerKval]
    end = first[-offsetLeft + 1 - spacerKval : -offsetLeft + 1]
    if (start in verificationIndexStartSpacer and end in verificationIndexEndSpacer and verificationIndexStartSpacer[start] == verificationIndexEndSpacer[end]):
        return True
    return False


# Reverse a given strand orientation
def ReverseOrientation(orientation):
    if (orientation == REGULAR):
        return REVCOMP
    else:
        return REGULAR
		

# Get the different types of orientation of a given set of reads
def GetDictStats(readOrientation, reads, readIndices, kmer):
    result = []
    for read in readIndices:
        state = readOrientation[reads[read]]
        if state != UNDEF and state not in result:
            result.append(state)
            if DEBUG and (len(result) > 1):
                print ("AMBIGUITY", reads[read], kmer)
    return result


# Orient a set of kmers and containing reads
def ReadKmerCanonization(reads, kmers, kmerToReads, kval):
    readOrientation = dict([(read, UNDEF) for read in reads])
    kmerOrientation = dict([(kmer, UNDEF) for kmer in kmerToReads])
    futureKmerOrientations = dict()
    bads = list()
    readlen = len(reads[0])
    carryon = True
    # Keep on orienting kmers and reads as long as there was an orientation in the prevoius round
    while (carryon):
        carryon = False
        for kmer in kmers:
            rc = RC(kmer)
            # Analyze kmers that were not already analyzed
            if (kmerOrientation[kmer] == UNDEF and kmerOrientation[rc] == UNDEF):
                # Get list of reads (for kmer and its reverse comlement)
                listOfReads = kmerToReads[kmer]
                listOfRCReads = kmerToReads[rc]
                # Verify that the no read contains the kmer and its reverse complement 
                commons = [read for read in listOfReads if read in listOfRCReads]
                if DEBUG and len(commons) > 0:
                    print ("BAD CANONIZATION: reads containing both kmer and RC:", kmer, rc)
                    continue
                # Get the orientations of reads containing the kmer (and its reverse complement)
                states = GetDictStats(readOrientation, reads, listOfReads, kmer)
                rcStates = GetDictStats(readOrientation, reads, listOfRCReads, rc)
                # Verify that no contradiction occurs in previous orientations of reads containig the kmer (or its reverse comlement)
                if (len(states) > 1 or len(rcStates) > 1):
                    if DEBUG:
                        print("BAD CANONIZATION: contradiction in orientation process")                
                        print (kmer, rc)
                    bads.append(kmer)
                    bads.append(rc)
                    continue
                # Verify that previous orientations of kmer and its reverse complement are not the same
                elif len(states) == 1 and len(rcStates) == 1 and states[0] == rcStates[0]:
                    if DEBUG:
                        print("BAD CANONIZATION: kmer and its reverse complement has the same orientation")
                        print (kmer, rc)
                    bads.append(kmer)
                    bads.append(rc)
                    continue
                # Orient kmer and reverse complement if possible
                if (len(states) == 0 and len(rcStates) == 0):
                    # Check if kmer can be oriented according to previously oriented kmers (actually, should not happen...)
                    if (kmer in futureKmerOrientations and rc in futureKmerOrientations):
                        # Verify no contradiction occurs in previous orientation of the kmer and its reverse complement
                        if (DEBUG and futureKmerOrientations[kmer] == futureKmerOrientations[rc]):
                            print("BAD CANONIZATION: contradiction in future kmer orientation")
                            continue
                        # Use previous orientation if possible
                        else:						
                            states.append(futureKmerOrientations[kmer])
                            rcStates.append(futureKmerOrientations[rc])                            						
                    # Determine orientation by current kmer (only for first root)
                    elif (len(futureKmerOrientations) == 0):	
                        states.append(REGULAR)
                        rcStates.append(REVCOMP)
                    # Orientation can not be determined by previously oriented kmers AND it is not the root oriented kmer - wait for others
                    else:
                        continue
                carryon = True
                # Make sure both kmer and its reverse complement are oriented 
                if (len(states) == 0):
                    states.append(ReverseOrientation(rcStates[0]))
                if (len(rcStates) == 0):
                    rcStates.append(ReverseOrientation(states[0]))
                # Orient kmer (mark kmer orientation, relevant reads orientation and  determine orientation of other kmers in those reads)
                kmerOrientation[kmer] = states[0]           
                for readIndex in listOfReads:
                    read = reads[readIndex]
                    readOrientation[read] = states[0]
                    for otherKmer in kmers:
                        if (readIndex in kmerToReads[otherKmer]):
                            if (otherKmer not in futureKmerOrientations):
                                futureKmerOrientations[otherKmer] = states[0]
                                futureKmerOrientations[RC(otherKmer)] = rcStates[0]
                # Orient reverse comlement kmer (in the same process)
                kmerOrientation[rc] = rcStates[0]
                for readIndex in listOfRCReads:
                    read = reads[readIndex]
                    readOrientation[read] = rcStates[0]
                    for otherKmer in kmers:
                        if (readIndex in kmerToReads[otherKmer]):
                            if (otherKmer not in futureKmerOrientations):
                                futureKmerOrientations[otherKmer] = rcStates[0]
                                futureKmerOrientations[RC(otherKmer)] = states[0]
    # Alert if orientation problems occured
    if (DEBUG and len(bads) > 0):
        print (len(kmerToReads))
        print (len(bads), bads)
        for kmer in kmers:
            print (kmer, kmerOrientation[kmer])
    return (bads, readOrientation, kmerOrientation)
	

# Obtain kmer data regarding a subset of the reads (that appear in a given dictionary)
# Construct a mapping between kmers and list of containing reads
def UpdateKmerData(reads, kval, stats):
    readindex = 0
    for read in reads:
        if (readindex % 100000 == 0):
            stdout.write("\r Updating k-mer data. processed: %d" % readindex)
            stdout.flush()
        for index in range(len(read) - kval + 1):
            kmer = read[index:index + kval]
            if (kmer in stats):
                occurences = stats[kmer]
                if (len(occurences) == 0 or occurences[-1] != readindex):
                    occurences.append(readindex)
                stats[kmer] = occurences
        readindex += 1
    stdout.write("\r Updating k-mer data. processed: %d" % readindex)
    stdout.flush()
    stdout.write("\n")
    return


# Cluster kmers according to intersection of the sets of containing reads
def ClusterKmers(kmersToReads):
    clusters = dict()
    # Initial clusters - one for every kmer
    identifier = 0
    for ((kmer, rckmer), (kmerReads, rcKmerReads)) in kmersToReads.items():
        kmerSet = set()
        kmerSet.add(kmer)
        kmerSet.add(rckmer)
        clusters[identifier] = [kmerSet, set(kmerReads).union(set(rcKmerReads))]
        identifier += 1
    merge = True
    # Keep merging until no two clusters can be merged
    while (merge):
        kmers = list(clusters.keys())
        merge = False
        firstClusterIdentifier = None
        secondClusterIdentifier = None
        for firstClusterIndex in range(len(kmers)):
            firstClusterIdentifier = kmers[firstClusterIndex]
            for secondClusterIndex in range(firstClusterIndex + 1, len(kmers)):
                secondClusterIdentifier = kmers[secondClusterIndex]
                if (len(clusters[firstClusterIdentifier][1].intersection(clusters[secondClusterIdentifier][1])) > MINIMUM_READS_CLUSTER_MERGE):
                    merge = True
                    break
            if (merge):
                break
        if (merge):
            clusters[firstClusterIdentifier][0] = clusters[firstClusterIdentifier][0].union(clusters[secondClusterIdentifier][0])
            clusters[firstClusterIdentifier][1] = clusters[firstClusterIdentifier][1].union(clusters[secondClusterIdentifier][1])
            del clusters[secondClusterIdentifier]
    return clusters


# Derive consensus sequence(s) for a repeat using a fast multiple alignment of all relevant reads (repeats)
def RepeatByMultipleAlignment(repeats, kmerShiftIndex, kmers, data, overlapsCache, minOverlap, readLen):
    if (len(repeats) == 0):
        return repeats
    # Multiple alignment of all repeats
    aligner = SimpleMultiple2([repData[0] for repData in repeats],KVAL, readLen, kmerShiftIndex)
    aligner.align()
    matrix = aligner.alignment
    if (DEBUG):
        for algn in matrix:
            print (algn)
    alignmentLength = aligner.maxLength
    repeatNumber = len(matrix)
    # Derive the multiple alignment statistics
    stats = dict()
    for baseIndex in range(alignmentLength):
        stats[baseIndex] = dict()
        for repIndex in range(repeatNumber):
            currBase = matrix[repIndex][baseIndex]
            stats[baseIndex][currBase] = stats[baseIndex].get(currBase, 0) + repeats[repIndex][1]
    if (DEBUG):
        print (stats)
    # Derive consensus sequences(s) for the repeat from multiple alignment statistics
    repeat = []
    dominants = []
    total = repeatNumber
    inRepeatMode = False
    for baseIndex in range(alignmentLength):
        # Remove spaces
        currTotal = total - stats[baseIndex].get(SPACE, 0)
        if (SPACE in stats[baseIndex]):
            stats[baseIndex].pop(SPACE)
        # Check for dominant bases in that position
        bases = []
        # A case of one dominant repeat base 
        dominant = False
        multipleBase = False
        if (DEBUG):
            print ("***", baseIndex, total, currTotal, stats[baseIndex])
        for base in stats[baseIndex]: 
            if (stats[baseIndex][base] > currTotal * CONSENSUS_FRACTION_DOMINANT):
                bases = [base]
                dominant = True
                if (DEBUG):
                    print ("Dominant", stats[baseIndex], base)
                break
        dominants.append(dominant)
        # A case of several (two) repeat base options
        if (not dominant):
            bases_sorted = sorted(stats[baseIndex], key=stats[baseIndex].get, reverse = True)
            if (stats[baseIndex][bases_sorted[0]] > currTotal * CONSENSUS_FRACTION_MULTIPLE_FIRST):
                if (stats[baseIndex][bases_sorted[1]] > currTotal * CONSENSUS_FRACTION_MULTIPLE_SECOND):
                        if ((len(stats[baseIndex])) < 3 or (stats[baseIndex][bases_sorted[2]] < currTotal * CONSEN_FRAC_MULT_MAX_THIRD)):
                            bases = [bases_sorted[0], bases_sorted[1]]
                            multipleBase = True
                            if (DEBUG):
                                print ("Multiple", stats[baseIndex], bases)
        if (DEBUG and not dominant and not multipleBase):
            print ("None", stats[baseIndex])
        # If there are dominant base(s) - add them
        if (len(bases) > 0):
            currRepeatBases = []
            for base in bases:
                currRepeatBases.append(base)
            repeat.append(currRepeatBases)
        else:
            repeat.append([SPACE])
    # Find longest sub sequence in the derived repeat
    startDR = currStart = endDR = currEnd = -1
    longestSequenceSoFar = currSequenceLength = 0
    isInSeq = False
    for index in range(len(repeat)):
        if (SPACE not in repeat[index]):
            if (not isInSeq):
                currStart = currEnd = index
                isInSeq = True
            else:
                currEnd = index
        else:
            if (isInSeq):
                isInSeq = False
                currSequenceLength = currEnd - currStart + 1
                if (currSequenceLength > longestSequenceSoFar):
                    longestSequenceSoFar = currSequenceLength
                    startDR = currStart
                    endDR = currEnd
    repeat = repeat[startDR : endDR + 1]
    # Convert repeat to string
    if (len(repeat) < MINIMUM_REPEAT_LENGTH or len(repeat) > MAXIMUM_REPEAT_LENGTH):
        return None
    repeatString = ""
    numberNonDominants = 0
    for repBase in repeat:
        if (len(repBase) == 1):
            repeatString += repBase[0]
        else:
            repeatString += "{"
            for option in repBase:
                repeatString += option + " "
            repeatString += "}"
            numberNonDominants += 1
    if (numberNonDominants > 9 or numberNonDominants > len(repeat) / 3):
        return None
    return repeatString
                   
 
# Determine a relative shift for every kmer according to a given root read
def ConstructShiftIndex(root, kmers):
    kmerShiftIndex = dict()
    for kmer in kmers:
        # kmer appears as it is in one of the initial DRs
        result = root.find(kmer)
        if (result > -1):
            kmerShiftIndex[kmer] = result       
            continue
        # kmer appeas only from the second base (flanking first base)
        result = root.find(kmer[1:])
        if (result > 0):
            kmerShiftIndex[kmer] = result - 1
            continue
        # kmer appeas only up to the last base (flanking last base)
        result = root.find(kmer[:-1])
        if (result > -1 and result <= len(root) - len(kmer)):
            kmerShiftIndex[kmer] = result
            continue
        # kmer does not appears exactly in one of the DRs
        aligner = LocalAligner(kmer, root)
        alignment = aligner.align()
        if (alignment.score > KVAL - 4):
            kmerShiftIndex[kmer] = alignment.shiftFirst
        else:
            #kmerShiftIndex[kmer] = 100
            pass
    return kmerShiftIndex
	
logger = open("log-crispr.txt", "a")
logger.write("\n" + "========================" + "\n")
logger.write("crispr Finder\n")
logger.write(FILE + " " + str(KVAL) + " " + str(threshold) + "\n")
startTime = time()

# Get all reads using a generator (memory efficient) 
print "Getting reads..."
reads1 = time()
readLen = GetReadsLength(DIRECTORY, FILE)
reads = readFileReadsThin(DIRECTORY + FILE)
reads2 = time()

# Get kmer frequency statistics: kmer => occurences
# Currently only kmers that appear once in a read are taken.
stats1 = time()
stats = GetKmerStats(reads, KVAL)
stats2 = time()
print "Calculating k-mer statistics took", stats2 - stats1, "seconds"

# Refine statistics to consider only frequent kmers
# Construct a mapping kmer => list of indexes of frequent kmers (and rverse complements as well).
# The number of occurences of a kmer and its reverse comlement are taken together.
# All lists are empty after invoking the method.
print "Taking frequent kmers..."
refinestats1 = time()
refinedStats = GetRefineStats(stats, threshold)
refinestats2 = time()
print "Taking frequent k-mers statistics took", refinestats2 - refinestats1, "seconds"
print "Number of frequent kmers is ", len(refinedStats)
logger.write("Number of frequent kmers is " + str(len(refinedStats)) + " " + str(len(refinedStats)/2) + "\n")
del stats

# Obtain the indices of the relevant reads (containing frequent kmers)
print "Taking relevant reads..."
refinereads1 = time()
refinedReads = GetRelevantReads(DIRECTORY, FILE, KVAL, refinedStats)
refinereads2 = time()
print "Taking relevant reads took", refinereads2 - refinereads1, "seconds"
print "Number of relevant reads is ", len(refinedReads)
logger.write ("Number of relevant reads is " + str(len(refinedReads)) + "\n")
del reads

# Calculate minimum overlap (consider reads length)
minOverlap = int(readLen * MIN_OVERLAP_FRACTION)
# Store the number of original relevant reads 
relevantReadsNumber = len(refinedReads)

# Get kmer data: every frequent kmer and the reads it appears in
# Update the data structure of kmer => list of indexes (that is currently with empty lists).
print "Update kmer data and duplicated reads"
updatekmers1 = time()
UpdateKmerData(refinedReads, KVAL, refinedStats)
updatekmers2 = time()

# Duplicate reads with their reverse complements
# We then have a list of reads r(1) .. r(n) and immediatley after we have RC(r(1)) .. RC(r(n))
duplicate1 = time()
refinedReads = DuplicateReads(refinedReads)
duplicate2 = time()
print "Updating data took", duplicate2 - updatekmers1, "seconds"

# Create index of hash values of prefixes for overlap queries (for reads and their rverse complementss)
hash1 = time()
data = ComputeHashValues(refinedReads, minOverlap, MOD)
hash2 = time()
print "Computing overlap index took", hash2 - hash1, "seconds"

# Cach of {kmer => {prefix length => overlapping indexes}}
overlapsCache = dict()
crisprKmers = []
print "Analyzing kmers..."
kindex = 0
startTimeKA = time()

# Match kmer pairs (kmer and its reversce complement)
match1 = time()
pairs = MatchKmers(refinedStats)
match2 = time()
gapt = 0
bapt = 0

# Work on every pair of frequent kmer and its reverse complement
for (kmer, rckmer) in pairs.items():
    # Get indices of reads containing thekmer and of reads containig the reverse complement
    readIndices = refinedStats[kmer]
    rcReadIndices = []
    if (rckmer != None):
	    rcReadIndices = refinedStats[rckmer]
    if (kindex % 1000 == 0):
        stdout.write("\r Analyzing kmers. processed: %d" % kindex)
        stdout.flush()		
    # Aggregate all indices: those of the kmer and those of the reverse complement 
    # (shifted towards the indices of the reverse complement reads)
    # This is done to have all reads from the same strand. 
    # Since this is a single kmer analysis, the kmer itself is chose as the right strand
    allIndices = []
    for regularIndex in readIndices:
        allIndices.append(regularIndex)
    for rcIndex in rcReadIndices:
        allIndices.append(relevantReadsNumber + rcIndex)
		
    # Analyze the kmer according to all reads and indices of kmer and revere complement (adjusted to reverse complement reads)
    apt1 = time()
    spacers = AnalyzeKmer(kmer, allIndices, refinedReads, data, overlapsCache, minOverlap)
    apt2 = time()
    kindex += 1
	
    # A frequent kmer is added to a list of tupples: (frequen kmer, frequent rverse complement)
    if (len(spacers) > 1):
        #print (kmer, spacers)
        crisprKmers.append((kmer, rckmer))
        gapt += apt2 - apt1
    else:
        bapt += apt2 - apt1
endTimeKA = time()
stdout.write("\r Getting relevant reads. processed: %d" % kindex)
stdout.flush()
stdout.write("\n")
print "Analyzing k-mers took", endTimeKA - startTimeKA, "seconds"

# Patition kmers into clusters
print "Clustring kmers..."
# Clustering should be done according to all reads that contain a kmer.
# This is due to the fact that the maximl paths of two consecutive kmers can potentially differ.
# Update the frequent kmers data with the lists of containing reads
crisprKmersData = dict([((kmerData[0], kmerData[1]), (refinedStats[kmerData[0]], refinedStats[kmerData[1]])) for kmerData in crisprKmers])
# Cluster kmers according to containing reads
cluster1 = time()
clusters = ClusterKmers(crisprKmersData)
cluster2 = time()
canon = 0;
distil = 0

for (identifier, (kmers, allReads)) in clusters.items():
    # Find orientation of all kmers and reads in a cluster
    canon1 = time()
    (bads, readOrientations, kmerOrientations) = ReadKmerCanonization(refinedReads, kmers, refinedStats, KVAL)
    canon2 = time()
    canon += canon2 - canon1
    if DEBUG and len(bads) > 0:
        print ("There are orientation problems.")
    # Orient the kmers according to the orientation algorithm
    orientedKmers = []
    for kmer in kmers:
        if (kmerOrientations[kmer] == REVCOMP):
            orientedKmers.append(RC(kmer))
        else:
            orientedKmers.append(kmer)
    # Orient the reads according to the orientation algorithm
    clusterReads = []
    for readIndex in allReads:
        if (readOrientations[refinedReads[readIndex]] == REVCOMP):
            clusterReads.append(RC(refinedReads[readIndex]))
        else:
            clusterReads.append(refinedReads[readIndex])
    # Construct list of (read, multiplicity) where multiplicity always equals one (no longer used)
    clusterReads = [[read,1] for read in clusterReads]
    # Choose the seed read to align accordingly (the read that contains the most kers)
    seedIndex = 0
    seedScore = -1
    for readIndex in range(len(clusterReads)):
        possibleSeed = clusterReads[readIndex][0]
        score = len([kmer for kmer in orientedKmers if kmer in possibleSeed])
        if (score > seedScore):
            seedScore = score
            seedIndex = readIndex
    # Switch seed read with first read
    temp = clusterReads[0]
    clusterReads[0] = clusterReads[seedIndex]
    clusterReads[seedIndex] = temp
    if (DEBUG):
        print ("SEED: ", clusterReads[seedIndex])
	
    # Construct shift index: for every kmer determine its relative location in the initial DR (seed read)
    kmerShiftIndex = ConstructShiftIndex(clusterReads[0][0], orientedKmers)
    if (DEBUG):
        print (clusterReads[0][0])
        for k in kmerShiftIndex:
            print (k, kmerShiftIndex[k])
    # Align all reads by kmer shift from inital repeats
    # in every read we search for the first kmer that also appears in the seed and align according to it
    distil1 = time()
    finalRepeat = RepeatByMultipleAlignment(clusterReads, kmerShiftIndex, orientedKmers, data, overlapsCache, minOverlap, readLen)
    distil2 = time()
    distil = distil2 - distil1
    if (finalRepeat != None):
        print "FINAL REPEAT:\n", finalRepeat, RC(finalRepeat)
        logger.write(finalRepeat + " " + RC(finalRepeat) + "\n")

endTime = time()
print "Total time:" , endTime - startTime, "seconds"

logger.write("STATS" + "\n")
logger.write("Reads: " + str(reads2 - reads1) + "\n")
logger.write("Stats: " + str(stats2 - stats1) + "\n")
logger.write("Refined Stats: " + str(refinestats2 - refinestats1) + "\n")
logger.write("Refined Reads: " + str(refinereads2 - refinereads1) + "\n")
logger.write("Update Kmers: " + str(updatekmers2 - updatekmers1) + "\n")
logger.write("Duplication: " + str(duplicate2 - duplicate1) + "\n")
logger.write("Hash Computation: " + str(hash2 - hash1) + "\n")
logger.write("Match Kmers: " + str(match2 - match1) + "\n")
logger.write("Good Analysis : " + str(gapt) + "\n")
logger.write("Bad Analysis: " + str(bapt) + "\n")
logger.write("Clustering: " + str(cluster2 - cluster1) + "\n")
logger.write("Canonization: " + str(canon) + "\n")
logger.write("Alignment: " + str(distil) + "\n")
logger.write("Total " + str(endTime - startTime) + "\n")
logger.close()

