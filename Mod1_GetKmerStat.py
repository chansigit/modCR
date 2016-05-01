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

# **********************************************************************************************************************
# *********************************************   Computation Section   ************************************************
# **********************************************************************************************************************
import json
import uuid
import logging
# 定义辅助参数
RESULTDIR= str(dictionaryArguments["-retdir"]) if "-retdir" in dictionaryArguments else "~/"
LOGDIR   = str(dictionaryArguments["-log"])    if "-log"    in dictionaryArguments else RESULTDIR
TASKNAME = str(dictionaryArguments["-tname"])  if "-tname"  in dictionaryArguments else "default"

# 生成日志文件名和输出文件名
uuidstr=str(uuid.uuid1())
KmerStats_Filename = "KmerStats_"+TASKNAME+"_"+uuidstr+".json"
KmerStats_Filename = os.path.join(RESULTDIR, KmerStats_Filename)
Log_Filename       = "KmerStats_"+TASKNAME+"_"+uuidstr+".log"
Log_Filename       = os.path.join(LOGDIR, Log_Filename)

# 写入任务信息
logging.basicConfig(level = logging.DEBUG, datefmt = '%a, %d %b %Y %H:%M:%S', filename = Log_Filename, filemode = 'w',
                    format = '%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
logtitle = "Module GetKmerStat - "+TASKNAME
logging.info(logtitle)
logging.info("Parameter K=" + str(KVAL) + ", T=" +str(threshold))
logging.info("Logging in "+str(Log_Filename))
logging.info("Results in "+str(KmerStats_Filename))

# 导入数据
logging.info("File Loading begins")
t1_load = time()
readLen = GetReadsLength(DIRECTORY, FILE)
reads = readFileReadsThin(DIRECTORY + FILE)
t2_load = time()
logging.info("File Loading Finished, taking " + str(t2_load-t1_load) + " seconds")

# 计算
logging.info("Getting Kmer Statistics begins")
t1_kmerStat = time()
stats = GetKmerStats(reads, KVAL)
t2_kmerStat = time()
logging.info("Getting Kmer Statistics ends, taking "+ str(t2_kmerStat-t1_kmerStat) +" seconds")

# 持久化结果
logging.info("Generating JSON file begins")
fp = open(KmerStats_Filename,"w")
t1_json =time()
json.dump(stats, fp)
t2_json =time()
fp.close()
logging.info("Generating JSON file ends, taking "+ str(t2_json-t1_json) +" seconds")