# This class represents an entry in the alignment matrix
# Entry consists of a score and an operation by which this entry was achieved.
from Aligners.Global import *;

class Entry:

    def __init__(self, score = DEFAULT_SCORE , operation = NO_OPERATOR):
        self.similarity = score;
        self.op = operation;
        
    def score(self):
        return self.similarity;

    def operation(self):
        return self.op;

    def setOperation(self, operation):
        self.op = operation;
        
    def setScore(self, score):
        self.similarity = score;    

    def final(self):
        return (self.op == NO_OPERATOR);
    
    def __lt__(self, other):
        return self.similarity < other.similarity;

    def __le__(self, other):
        return self.similarity <= other.similarity;

    def __eq__(self, other):
        return self.similarity == other.similarity;

    def __ne__(self, other):
        return self.similarity != other.similarity;

    def __gt__(self, other):
        return self.similarity > other.similarity;

    def __ge__(self, other):
        return self.similarity >= other.similarity;
    
    def __repr__(self):
        string = str(self.similarity) + " " + str(self.op);
        return string;    
