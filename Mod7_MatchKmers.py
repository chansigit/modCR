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
# Reverse complement a sequence
def RC(seq):
    return ''.join([RC_DICT[base] for base in seq[::-1]])

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
# **********************************************************************************************************************
# *********************************************   Computation Section   ************************************************
# **********************************************************************************************************************
# python Mod7_MatchKmers.py c:\data\seq.fasta -k 20 -t 40 -retdir c:\data -log c:\data\ -tname may1test -stat C:\data\Mod4_RefinedStats_may1test_900dbb21-11c4-11e6-8c9c-ec55f98094e4.json

import cPickle
import uuid
import logging
# 定义辅助参数
RESULTDIR= str(dictionaryArguments["-retdir"]) if "-retdir" in dictionaryArguments else "~/"
LOGDIR   = str(dictionaryArguments["-log"])    if "-log"    in dictionaryArguments else RESULTDIR
TASKNAME = str(dictionaryArguments["-tname"])  if "-tname"  in dictionaryArguments else "default"
refinedStatsFile=str(dictionaryArguments["-stat"])  if "-stat"  in dictionaryArguments else "refinedstat.pkl"

# 生成日志文件名和输出文件名
uuidstr=str(uuid.uuid1())
pair_Filename = "Mod7_MatchedPair_"+TASKNAME+"_"+uuidstr+".pkl"  #
pair_Filename = os.path.join(RESULTDIR, pair_Filename)            #
Log_Filename  = "Mod7_MatchedPair_"+TASKNAME+"_"+uuidstr+".log"
Log_Filename  = os.path.join(LOGDIR, Log_Filename)

# 写入任务信息
logging.basicConfig(level = logging.DEBUG, datefmt = '%a, %d %b %Y %H:%M:%S', filename = Log_Filename, filemode = 'w',
                    format = '%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
logtitle = "Module RefinedStat - "+TASKNAME
logging.info(logtitle)
logging.info("Parameter K=" + str(KVAL) + ", T=" +str(threshold))
logging.info("Logging in "+str(Log_Filename))
logging.info("Results in "+str(pair_Filename))

# 导入数据
logging.info("File Loading begins")
t1_load = time()
refinedStats = cPickle.load(open(refinedStatsFile, 'r'))
t2_load = time()
logging.info("File Loading Finished, taking " + str(t2_load-t1_load) + " seconds")

# 计算
logging.info("Matching Kmer Pairs begins")
t1_pairMatch = time()
pairs = MatchKmers(refinedStats)
t2_pairMatch = time()
logging.info("Matching Kmer Pairs ends, taking "+ str(t2_pairMatch-t1_pairMatch) +" seconds")

# 持久化结果
logging.info("Serialization Begins")
fp = open(pair_Filename,"w")
t1_serial =time()
cPickle.dump(pairs, fp)
t2_serial =time()
fp.close()
logging.info("Serialization Ends, taking "+ str(t2_serial-t1_serial) +" seconds")