# Pairwise alignment of strings
from Aligners.Aligner import *;
from Aligners.Global import *;
from Aligners.Entry import *;
from Aligners.Alignment import *;

class Pairwise(Aligner):
    
    # Base condition is as the prefix size (either row or col)
    def baseCondition(self, row, col):
        return Entry(max(row, col), NO_OPERATOR);
    
    # Get coordinates of cell in the last column and row
    def maxEntryCoordinates(self):
        return (len(self.first), len(self.second))
    
    def suffixAlignment(self, maxEntryCoordinates):
        return ([], [])
    
    def prefixAlignment(self, beginCoordinates):
        prefixFirst = []
        prefixSecond = []
        if (beginCoordinates[0] > 0):
            for index in range(beginCoordinates[0]):
                prefixFirst.append(self.first[index])
                prefixSecond.append(SPACE)
        elif (beginCoordinates[1] > 0):
            for index in range(beginCoordinates[1]):
                prefixSecond.append(self.second[index])
                prefixFirst.append(SPACE)
        return (prefixFirst, prefixSecond)
    
    
