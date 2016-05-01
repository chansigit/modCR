# Pairwise alignment of overlapping strings
from Aligners.Aligner import *;
from Aligners.Global import *;
from Aligners.Entry import *;
from Aligners.Alignment import *;

class EndsFree(Aligner):
    
    # Base condition is always zero for ends free alignment
    def baseCondition(self, row, col):
        return Entry(0, NO_OPERATOR);
    
    
    # Get coordinates of entry with the highest score in last column/row
    def maxEntryCoordinates(self):
        # Initial values
        maxValue = float("-inf");
        maxCoordinates = (-1,-1);
        # Get the the last row
        lastRow = self.matrix[len(self.first)];        
        # Get the the last column
        lastCol = [row[len(self.second)] for row in self.matrix];
        rowIndex = len(self.first);
        # Go over the last row and maintain the coordinates of the maximum entry
        for colIndex in range(len(self.second) + 1):
            # Update data if a bigger value is encountered in the row
            if (lastRow[colIndex].score() > maxValue):
                maxValue = lastRow[colIndex].score();
                maxCoordinates = (rowIndex, colIndex);        
        # Go over the last column and maintain the coordinates of the maximum entry
        colIndex = len(self.second);
        for rowIndex in range(len(self.first) + 1):
            # Update data if a bigger value is encountered in the column
            if (lastCol[rowIndex].score() > maxValue):
                maxValue = lastCol[rowIndex].score();
                maxCoordinates = (rowIndex, colIndex);   
        # Return the coordinates of the maximum entry
        return maxCoordinates;    
    
    # Get the suffixes of the aligned strings
    def suffixAlignment(self, maxEntryCoordinates):
        # Calculate lengths
        sizeFirst = len(self.first);
        sizeSecond = len(self.second);

        # In case strings end together
        suffix1 = [];
        suffix2 = [];
        
        # In case the first string ended before the second one
        diff = sizeSecond - maxEntryCoordinates[1];
        if(diff > 0):
            suffix1 = [SPACE for index in range(diff)];
            suffix2 = [self.second[index] for index in range(maxEntryCoordinates[1], sizeSecond)];
                                                        
        # In case the second string ended before the first one
        diff = sizeFirst - maxEntryCoordinates[0];
        if(diff > 0):
            suffix1 = [self.first[index] for index in range(maxEntryCoordinates[0], sizeFirst)];
            suffix2 = [SPACE for index in range(diff)];

        # return the suffixes
        return (suffix1, suffix2);
    
  
    # Get the prefixes of the aligned strings
    def prefixAlignment(self, beginCoordinates):
        
        # extract coordinates
        counter1 = beginCoordinates[0];
        counter2 = beginCoordinates[1];
        
        # In case strings start together
        prefix1 = [];
        prefix2 = [];
        
        # In case the first string starts before the second one
        if(counter1 > 0):
            prefix1 = [self.first[index] for index in range(counter1)];
            prefix2 = [SPACE for index in range(counter1)];
        
        # In case the second string starts before the first one
        if(counter2 > 0):
            prefix1 = [SPACE for index in range(counter2)];
            prefix2 = [self.second[index] for index in range(counter2)];        

        # return the prefixes
        return (prefix1, prefix2);
    
#str1 = "abcdefghijklmnopq";
#str2 = "fhqjklq";
#c = EndsFree(str1, str2);
#print c.align()
