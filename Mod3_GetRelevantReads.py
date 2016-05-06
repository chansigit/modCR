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

# Get the length of the reads (assuming all have the same length)
def GetReadsLength(directory, fileName):
    handler = open(directory + fileName, "r")
    line = handler.readline().strip()
    while line.startswith(">"):
        line = handler.readline().strip()
    return len(line);

# **********************************************************************************************************************
# *********************************************   Computation Section   ************************************************
# **********************************************************************************************************************
# python Mod3_GetRelevantReads.py c:\data\seq.fasta -k 20 -t 40 -retdir c:\data -log c:\data\ -tname may1test -rfst C:\data\Mod2_RefinedStats_may1test_a4621a8f-11be-11e6-be68-ec55f98094e4.json

import cPickle
import uuid
import logging
# 定义辅助参数
RESULTDIR= str(dictionaryArguments["-retdir"]) if "-retdir" in dictionaryArguments else "~/"
LOGDIR   = str(dictionaryArguments["-log"])    if "-log"    in dictionaryArguments else RESULTDIR
TASKNAME = str(dictionaryArguments["-tname"])  if "-tname"  in dictionaryArguments else "default"
refineStatsFile=str(dictionaryArguments["-rfst"])  if "-rfst"  in dictionaryArguments else "refinestat.pkl"  #

# 生成日志文件名和输出文件名
uuidstr=str(uuid.uuid1())
RefinedReads_Filename = "Mod3_RefinedReads_" +TASKNAME+"_"+uuidstr+".pkl"  #
RefinedReads_Filename = os.path.join(RESULTDIR, RefinedReads_Filename)      #

readLen_Filename      = "Mod3_readLen_"       +TASKNAME+"_"+uuidstr+".pkl"  #
readLen_Filename      = os.path.join(RESULTDIR, readLen_Filename     )      #

minOverlap_Filename   = "Mod3_minOverlap_"    +TASKNAME+"_"+uuidstr+".pkl"  #
minOverlap_Filename   = os.path.join(RESULTDIR, minOverlap_Filename  )      #

relevantReadsNumber_Filename = "Mod3_relevantReadsNumber_" +TASKNAME+"_"+uuidstr+".pkl"  #
relevantReadsNumber_Filename = os.path.join(RESULTDIR, relevantReadsNumber_Filename)     #

Log_Filename          = "Mod3_RefinedReads_"+TASKNAME+"_"+uuidstr+".log"   #
Log_Filename          = os.path.join(LOGDIR, Log_Filename)

# 写入任务信息
logging.basicConfig(level = logging.DEBUG, datefmt = '%a, %d %b %Y %H:%M:%S', filename = Log_Filename, filemode = 'w',
                    format = '%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
logtitle = "Module RefinedReads - "+TASKNAME             #
logging.info(logtitle)
logging.info("Parameter K=" + str(KVAL) + ", T=" +str(threshold))
logging.info("Logging in "+str(Log_Filename))
logging.info("Results in "+str(RefinedReads_Filename)) #
logging.info("Results in "+str(readLen_Filename))      #
logging.info("Results in "+str(minOverlap_Filename))   #
logging.info("Results in "+str(relevantReadsNumber_Filename))   #

# 导入数据
logging.info("File Loading begins")
t1_load     = time()
refineStats = cPickle.load(open(refineStatsFile, 'r'))    #
t2_load     = time()
logging.info("File Loading Finished, taking " + str(t2_load-t1_load) + " seconds")

# 计算
logging.info("Refining Reads begins")             #
t1_refineReads = time()    #
refinedReads   = GetRelevantReads(DIRECTORY, FILE, KVAL, refineStats)     #
readLen        = GetReadsLength(DIRECTORY, FILE)      #
minOverlap     = int(readLen * MIN_OVERLAP_FRACTION)  #
relevantReadsNumber = len(refinedReads)           #
t2_refineReads = time()    #
logging.info("Refining Reads ends, taking "+ str(t2_refineReads-t1_refineReads) +" seconds")    #

# 持久化结果
logging.info("Serialization begins")
t1_serial =time()

fp = open(RefinedReads_Filename,"w")   #
cPickle.dump(refinedReads, fp)            #
fp.close()

fp = open(readLen_Filename,"w")   #
cPickle.dump(readLen, fp)            #
fp.close()

fp = open(minOverlap_Filename,"w")#
cPickle.dump(minOverlap, fp)            #
fp.close()

fp = open(relevantReadsNumber_Filename,"w")   #
cPickle.dump(relevantReadsNumber, fp)            #
fp.close()

t2_serial =time()
logging.info("Serialization ends, taking "+ str(t2_serial-t1_serial) +" seconds")