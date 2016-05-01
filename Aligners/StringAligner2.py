from Aligners.Alignment import *;
from Aligners.Global import *;

class StringAligner2():

    # Constructor
    def __init__(self, first, second, kval, readlen, kmerShiftIndex):
        # Store relevant strings
        self.first = first;
        self.second = second;
        self.kval = kval;
        self.readlen = readlen
        self.kmerShiftIndex = kmerShiftIndex

    def align(self):
        kmerShiftInSecond = -1
        for index in range(self.readlen - self.kval + 1):
            kmer = self.second[index : index + self.kval]
            if (kmer in self.kmerShiftIndex):
                kmerShiftInFirst = self.kmerShiftIndex[kmer]
                kmerShiftInSecond = index
                break
        if (kmerShiftInSecond < 0):
            print ("Error in alignment!!")
            return None
        diff = kmerShiftInFirst - kmerShiftInSecond
        if (diff > 0):
            aFirst = []
            for base in self.first:
                aFirst.append(base)
            for i in range(diff):
                aFirst.append("-")
            aSecond = []
            for i in range(diff):
                aSecond.append("-")
            for base in self.second:
                aSecond.append(base)
            align = Alignment(self.readlen - diff, (aFirst , aSecond) , 0);
        elif (diff < 0):
            diff = -diff
            aSecond = []
            for base in self.second:
                aSecond.append(base)
            for i in range(diff):
                aSecond.append("-")
            aFirst = []
            for i in range(diff):
                aFirst.append("-")
            for base in self.first:
                aFirst.append(base)
            align = Alignment(self.readlen - diff, (aFirst , aSecond) , diff);
        else:
            aFirst = []
            for base in self.first:
                aFirst.append(base)
            aSecond = []
            for base in self.second:
                aSecond.append(base)
            align = Alignment(self.readlen, (aFirst , aSecond) , 0);
        return align        
    

    def matches(self, first, second):
        return len([index for index in range(len(first)) if first[index] == second[index]])

    def align2(self):       
        found1 = False
        found2 = False
        for shift in range(1, self.readlen - self.kval + 1):
            if (self.matches(self.first[shift:], self.second[:self.readlen - shift]) > self.kval - 1):
                found1 = True
                break
            if (self.matches(self.second[shift:], self.first[:self.readlen - shift]) > self.kval - 1):
                found2 = True
                break
        if (found1):
             aFirst = []
             for base in self.first:
                 aFirst.append(base)
             for i in range(shift):
                 aFirst.append("-")
             aSecond = []
             for i in range(shift):
                 aSecond.append("-")
             for base in self.second:
                 aSecond.append(base)
             align = Alignment(self.readlen - shift, (aFirst , aSecond) , 0);
        elif (found2):
             aSecond = []
             for base in self.second:
                 aSecond.append(base)
             for i in range(shift):
                 aSecond.append("-")
             aFirst = []
             for i in range(shift):
                 aFirst.append("-")
             for base in self.first:
                 aFirst.append(base)
             align = Alignment(self.readlen - shift, (aFirst , aSecond) , shift);
        else:
             align = Alignment(0, ([] , []) , 0);
        return align

            



        

