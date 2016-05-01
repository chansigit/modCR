from abc import *;
# Inteface for pairwise alignment calculator
class IAligner:
    
    __metaclass__ = ABCMeta;
    
    @abstractmethod
    # Reset internal data
    def reset(self):
        return NotImplemented;   
    
    @abstractmethod
    # Get an alignment object between of two strings
    def align(self):
        return NotImplemented;
    
    @abstractmethod
    # Create a dynamic programming matrix for pairwise alignment
    def createMatrix(self):
        return NotImplemented;    
    
    @abstractmethod
    # Calculate the value of an entry
    def entryValue(self, row, col):
        return NotImplemented;
    
    @abstractmethod
    # get a base condition value
    def baseCondition(self, row, col):
        return NotImplemented;
    
    @abstractmethod
    # get a value of an inner entry in the matrix
    def reoccurenceValue(self, row, col):
        return NotImplemented;
    
    @abstractmethod
    # Get the coordinates of the entry that contains the final score
    def maxEntryCoordinates(self):
        return NotImplemented;
    
    @abstractmethod    
    # Get aligned strings starting from the final entry
    def traceback(self, finalEntryCoordinates):    
        return NotImplemented;
    
    @abstractmethod  
    # Get the suffixes of the aligned strings (as tuple)
    def suffixAlignment(self, maxEntryCoordinates):
        return NotImplemented;
    
    @abstractmethod  
    # Get the prefixes of the aligned strings (as tuple)
    def prefixAlignment(self, beginCoordinates):
        return NotImplemented;     
    
    @abstractmethod    
    # Compare two characters
    def compare(self, char1, char2):    
        return NotImplemented;    
    
    @classmethod
    def __subclasshook__(cls, C):
        return NotImplemented;