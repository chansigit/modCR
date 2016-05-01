# Operations
SUBSTITUTION_OPERATOR = 'S';
DELETION_OPERATOR = 'D';
INSERTION_OPERATOR = 'I';
EQUALITY_OPERATOR = "E";
NO_OPERATOR = "N";

# Character Constants
SPACE = "-";
EQUALITY = "=";
MISMATCH_PREFIX = "(";
MISMATCH_SUFFIX = ")";
JOKER = "X";

# Score
MATCH_SCORE = 3;
MISMATCH_SCORE = -1;
DEFAULT_SCORE = -1;
BAD_ALLIGNMENT_SCORE = float("-inf");
OPTIMAL_K = 13;
ERROR_RATE = 0.1
#THRESHOLD_FACTOR = 1.5;
THRESHOLD_FACTOR = 2;

# File Names
GENOME_FILE = "Genomes.txt";
PROMOTORS_FILE = "Promotors.txt";

# Regular Expressions
SEPERATOR_REGEXP = "=+";

# Get the minimum score of an alignment
def getMinScore(first, second):
    return 0;

# Read entire file content to a list
def readFileContent(filename):
    # Create an handler to the file
    fileHandler = open(filename, 'r');    
    data = [];
    # read lines one by one
    line = fileHandler.readline();
    while(line != ''):
        # trim the line
        line = line.strip();
        # Add the line to the content
        data.append(line);
        # Read the next line
        line = fileHandler.readline();
    # Close handler
    fileHandler.close();
    return data;