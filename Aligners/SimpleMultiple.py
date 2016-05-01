from Aligners.StringAligner import *;
from Aligners.Alignment import *;
from Aligners.Global import *;

class SimpleMultiple():

    def __init__(self, sequences, kval, readlen, kmerShiftIndex, initalRepeatShift):
        self.sequences = sequences
        self.seedShift = 0
        self.maxLength = 0
        self.number = len(sequences)
        self.alignment = []
        self.kval = kval
        self.readlen = readlen
        self.kmerShiftIndex = kmerShiftIndex
        self.initalRepeatShift = initalRepeatShift

    def align(self):
        if (self.number == 0):
            return 
        # Seed choice
        seed = self.sequences[0]
        if (self.number == 1):
            seedAsList = []
            for base in seed:
                seedAsList.append(base)                
            self.maxLength = len(seed)
            self.alignment.append(seedAsList)
            return 
        # Perform initial local alignment between a seed and another sequence
        aligner = StringAligner(seed, self.sequences[1], self.kval, self.readlen, self.kmerShiftIndex, self.initalRepeatShift)
        result = aligner.align()
        self.alignment.append(result.first())
        self.alignment.append(result.second())
        self.seedShift = result.shiftFirst
        self.maxLength = len(result.first())
        # Perform sequential local alignments between all sequences and the seed
        for index in range(2, self.number):
            aligner = StringAligner(seed, self.sequences[index], self.kval, self.readlen, self.kmerShiftIndex, self.initalRepeatShift)
            result = aligner.align()
            # Insert the new aligned sequence
            self.alignment.append(result.second())
            # Update the prefix of the alignments (as required)
            addition = 0
            additionAll = 0
            if (result.shiftFirst == 0):
                addition = self.seedShift
                for spaceIndex in range(self.seedShift):
                    self.alignment[index].insert(0, SPACE)
            elif (result.shiftFirst <= self.seedShift):
                addition = self.seedShift - result.shiftFirst
                for spaceIndex in range(self.seedShift - result.shiftFirst):
                    self.alignment[index].insert(0, SPACE)
            else:
                additionAll = result.shiftFirst - self.seedShift
                for prev in range(index):
                    for spaceIndex in range(result.shiftFirst - self.seedShift):
                        self.alignment[prev].insert(0, SPACE)
                self.seedShift = result.shiftFirst
            # Insert ending spaces (as required)
            if (self.maxLength + additionAll < len(result.first()) + addition):
                for prev in range(index):
                    for spaceIndex in range(len(result.first()) + addition - (self.maxLength + additionAll)):
                        self.alignment[prev].append(SPACE)
                self.maxLength = len(result.first()) + addition
            else:
                #print (self.maxLength, additionAll, len(result.first()), addition)
                for spaceIndex in range(self.maxLength + additionAll - (len(result.first()) + addition)):
                    self.alignment[index].append(SPACE)
                self.maxLength = self.maxLength + additionAll
        for algn in self.alignment:
            for elemIndex in range(len(algn)):
                if algn[elemIndex].startswith(MISMATCH_PREFIX):
                    algn[elemIndex] = algn[elemIndex][1]
            #print (algn)
        #print (self.seedShift)
        #print (self.maxLength)
