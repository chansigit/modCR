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
# **********************************************************************************************************************
# *********************************************   Computation Section   ************************************************
# **********************************************************************************************************************
# python Mod4_UpdateKmerData.py c:\data\seq.fasta -k 20 -t 40 -retdir c:\data -log c:\data\ -tname may1test -rfst C:\data\Mod2_RefinedStats_may1test_a4621a8f-11be-11e6-be68-ec55f98094e4.json -rfrd C:\data\Mod3_RefinedReads_may1test_54a09d91-11d4-11e6-802c-ec55f98094e4.json
import json
import uuid
import logging
# 定义辅助参数
RESULTDIR= str(dictionaryArguments["-retdir"]) if "-retdir" in dictionaryArguments else "~/"
LOGDIR   = str(dictionaryArguments["-log"])    if "-log"    in dictionaryArguments else RESULTDIR
TASKNAME = str(dictionaryArguments["-tname"])  if "-tname"  in dictionaryArguments else "default"
refinedStatsFile=str(dictionaryArguments["-rfst"])  if "-rfst"  in dictionaryArguments else "refinestat.json"  #
refinedReadsFile=str(dictionaryArguments["-rfrd"])  if "-rfrd"  in dictionaryArguments else "refineread.json"  #

# 生成日志文件名和输出文件名
uuidstr=str(uuid.uuid1())
RefinedStats_Filename = "Mod4_RefinedStats_"+TASKNAME+"_"+uuidstr+".json"  #
RefinedStats_Filename = os.path.join(RESULTDIR, RefinedStats_Filename)     #
Log_Filename          = "Mod4_RefinedStats_"+TASKNAME+"_"+uuidstr+".log"   #
Log_Filename          = os.path.join(LOGDIR, Log_Filename)

# 写入任务信息
logging.basicConfig(level = logging.DEBUG, datefmt = '%a, %d %b %Y %H:%M:%S', filename = Log_Filename, filemode = 'w',
                    format = '%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
logtitle = "Module UpdateKmer - "+TASKNAME             #
logging.info(logtitle)
logging.info("Parameter K=" + str(KVAL) + ", T=" +str(threshold))
logging.info("Logging in "+str(Log_Filename))
logging.info("Results in "+str(RefinedStats_Filename))   #

# 导入数据
logging.info("File Loading begins")
t1_load     = time()
refinedStats = json.load(open(refinedStatsFile, 'r'))    #
refinedReads = json.load(open(refinedReadsFile, 'r'))    #
t2_load     = time()
logging.info("File Loading Finished, taking " + str(t2_load-t1_load) + " seconds")

# 计算
logging.info("Update Kmer begins")          #
t1_updateKmer = time()
UpdateKmerData(refinedReads, KVAL, refinedStats)
t2_updateKmer = time()
logging.info("Update Kmer ends, taking "+ str(t2_updateKmer-t1_updateKmer) +" seconds")    #

# 持久化结果
logging.info("Generating JSON file begins")
fp = open(RefinedStats_Filename,"w")   #
t1_json =time()
json.dump(refinedStats, fp)            #
t2_json =time()
fp.close()
logging.info("Generating JSON file ends, taking "+ str(t2_json-t1_json) +" seconds")