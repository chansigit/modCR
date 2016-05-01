# Base class containing methods of pairwise alignment
# The file contains general operations regarding all kinds of alignment.
# Implementation of specific alignment methods is postponed to derived classes.
from Aligners.IAligner import *;
from Aligners.Entry import *;
from Aligners.Alignment import *;
from Aligners.Global import *;

class FullLocalAligner(IAligner):
    
    # Constructor
    def __init__(self, first, second):
        # Store relevant strings
        self.first = first;
        self.second = second;
        # Set initial score to zero
        self.score = 0;
        # Construct an empty dynamic programming matrix
        self.matrix = [];
        self.shiftFirst = -1
        self.shiftSecond = -1
        
        
    # Reset internal data
    def reset(self):
        self.score = 0;
        self.matrix = [];
    
    
    # Get an alignment object between of two strings
    def align(self):
        
        # Create the dynamic programming matrix
        self.createMatrix();
    
        # Find the proper maximum entry
        maxEntryCoordinates = self.maxEntryCoordinates();
    
        # Get the allignment score
        score = self.matrix[maxEntryCoordinates[0]][maxEntryCoordinates[1]].score();
    
        # Get the threshild for a good alignment
        threshold = getMinScore(self.first, self.second);
        # Construct the allignments, if necessary (and the reads are not identical)
        allignments = ("", "");
        if(score >= threshold):
            (shiftFirst, allignments) = self.traceback(maxEntryCoordinates);
        else:
            score = BAD_ALLIGNMENT_SCORE;
    
        # construct an allignment object
        allign = Alignment(score, allignments, shiftFirst);
        return allign;
    
    
    # Create a dynamic programming matrix for pairwise alignment
    def createMatrix(self):
        # Allocate memory for the dynamic programming matrix
        rows = len(self.first) + 1;
        cols = len(self.second) + 1;
        self.matrix = [[Entry()]*cols for row in range(rows)];
        
        # Fill all matrix entries
        for row in range(rows):
            for col in range(cols):
                self.matrix[row][col] = self.entryValue(row, col);      
    
    
    # Calculate the value of an entry
    def entryValue(self, row, col):
        # Base condition case
        if ((row == 0) or (col == 0)):
            result = self.baseCondition(row, col);
        # Reocuurence case
        else:
            result = self.reoccurenceValue(row, col);
        return result;
    
    
    # get a base condition value
    def baseCondition(self, row, col):
        return Entry(0, NO_OPERATOR);
    
    
    # get a value of an inner entry in the matrix
    def reoccurenceValue(self, row, col):
        
        # Get he relevant characters
        first = self.first[row - 1];
        second = self.second[col - 1];
        
        # Compare the last two chars (equality or substitution)
        baseEqualityScore = self.compare(first, second)
        substitute = (self.matrix[row - 1][col - 1]).score() + baseEqualityScore;
        # Decide the right operation (equality or substitution)
        if (baseEqualityScore == MATCH_SCORE):
            operation = EQUALITY_OPERATOR;
        else:
            operation = SUBSTITUTION_OPERATOR;
        subs = Entry(substitute, operation);
    
        # Compare the first character to space (deletion)
        deletion = (self.matrix[row - 1][col]).score() + self.compare(first, SPACE);
        dele = Entry(deletion, DELETION_OPERATOR);
    
        # Compare the second character to space (insertion)
        insertion = (self.matrix[row][col - 1]).score() + self.compare(SPACE, second);
        inse = Entry(insertion, INSERTION_OPERATOR);
		
		# Local alignemnt has an option of no op in the table itself as well
        noop = Entry(0, NO_OPERATOR)
                        
        # Get the maximum similarity
        #value = max(subs, dele, inse);
        value = max(inse, dele, subs, noop)
    
        return value;
    
    
    # Get the coordinates of the entry that contains the final score
    def maxEntryCoordinates(self):
        rows = len(self.first) + 1;
        cols = len(self.second) + 1;
        maxCoordinates = (0,0)
        maxScore = 0
        for row in range(rows):
            for col in range(cols):
                if (self.matrix[row][col].score() > maxScore):
                    maxCoordinates = (row, col)
                    maxScore = self.matrix[row][col].score()
        return maxCoordinates
    
    # Get aligned strings starting from the final entry
    def traceback(self, maxEntryCoordinates):
        
        # Construct the aligned strings suffixes
        suffixes = self.suffixAlignment(maxEntryCoordinates);
        str1 = suffixes[0];
        str2 = suffixes[1];
        
        # Get the final coordinates data
        rowCounter = maxEntryCoordinates[0];  
        colCounter = maxEntryCoordinates[1];
        
        # Walk through the matrix from the maximum entry till a final entry
        currEntry = self.matrix[rowCounter][colCounter];
        while(not currEntry.final()):
    
            # In case of an equality - take the both equal characters
            if(currEntry.operation() == EQUALITY_OPERATOR):
                str1.insert(0, self.first[rowCounter - 1]);
                str2.insert(0, self.second[colCounter - 1]);
                rowCounter = rowCounter - 1;
                colCounter = colCounter - 1;
                
            # In case of a substitution - take both characters with mismatch markers
            elif(currEntry.operation() == SUBSTITUTION_OPERATOR):
                str1.insert(0, MISMATCH_PREFIX + self.first[rowCounter - 1] + MISMATCH_SUFFIX);
                str2.insert(0, MISMATCH_PREFIX + self.second[colCounter - 1] + MISMATCH_SUFFIX);
                rowCounter = rowCounter - 1;
                colCounter = colCounter - 1;
    
            # In case of a deletion - take first character and a SPACE from the second
            elif(currEntry.operation() == DELETION_OPERATOR):
                str1.insert(0, self.first[rowCounter - 1]);
                str2.insert(0, SPACE);
                rowCounter = rowCounter - 1;
    
            # In case of an insertion - take second character and a SPACE from the first
            elif(currEntry.operation() == INSERTION_OPERATOR):
                str1.insert(0, SPACE);
                str2.insert(0, self.second[colCounter - 1]);
                colCounter = colCounter - 1;
    
            # Error in operation
            else:
                print ("Error In Operation!");
    
            # Update the current entry
            currEntry = self.matrix[rowCounter][colCounter];
    
        # Construct the aligned strings prefixes
        beginCoordinates = (rowCounter, colCounter);
        prefixes = self.prefixAlignment(beginCoordinates);
    
        # Construct the final aligned strings
        allignedStr1 = prefixes[0];
        allignedStr1.extend(str1);
        allignedStr2 = prefixes[1];
        allignedStr2.extend(str2);

        self.shiftFirst = rowCounter
        self.shiftSecond = colCounter

        return (max(0, beginCoordinates[1] - beginCoordinates[0]), (allignedStr1, allignedStr2));
    

    # Get the suffixes of the aligned strings
    def suffixAlignment(self, maxEntryCoordinates):
        # Calculate lengths
        sizeFirst = len(self.first);
        sizeSecond = len(self.second);
        
        # In case the first string ended before the second one
        suffix1 = [self.first[index] for index in range(maxEntryCoordinates[0], sizeFirst)];
        suffix2 = [self.second[index] for index in range(maxEntryCoordinates[1], sizeSecond)];

        # Complete the aligment with spaces 
        diff1 = sizeFirst - maxEntryCoordinates[0];
        diff2 = sizeSecond - maxEntryCoordinates[1];
        for index in range(diff1, diff2):
            suffix1.append(SPACE)
        for index in range(diff2, diff1):
            suffix2.append(SPACE)

        # return the suffixes
        return (suffix1, suffix2);
    
  
    # Get the prefixes of the aligned strings
    def prefixAlignment(self, beginCoordinates):
        # Calculate lengths
        sizeFirst = len(self.first);
        sizeSecond = len(self.second);
        
        # In case the first string ended before the second one
        prefix1 = [self.first[index] for index in range(beginCoordinates[0])];
        prefix2 = [self.second[index] for index in range(beginCoordinates[1])];

        # Complete the aligment with spaces 
        for index in range(beginCoordinates[0], beginCoordinates[1]):
            prefix1.insert(0, SPACE)
        for index in range(beginCoordinates[1], beginCoordinates[0]):
            prefix2.insert(0, SPACE)

        # return the suffixes
        return (prefix1, prefix2); 
    
   
    # Compare two characters
    def compare(self, char1, char2):
        result = 0;
        if(char1 == char2 or ((char1 != SPACE and char2 != SPACE) and (char1 == JOKER or char2 == JOKER))):
            result = MATCH_SCORE;
        else:
            result = MISMATCH_SCORE;
        return result;
    
    # Representation operator
    def __repr__(self):
        result =  self.first + "\n";
        result =  result + self.second + "\n";
        result =  result + repr(self.matrix) + "\n";
        return result;
            
