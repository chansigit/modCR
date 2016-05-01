from Aligners.Global import *;

# This class represents an alignment (score and two strings aligned)
class Alignment:

    def __init__(self, score, data, shiftFirst = -1):
        self.score = score;
        self.firstString = data[0];
        self.secondString = data[1];
        self.shiftFirst = shiftFirst

    def similarity(self):
        return self.score;

    def first(self):
        return self.firstString;

    def second(self):
        return self.secondString;
    
    def __repr__(self):
        string = str(self.score) + "\n";
        # Construct the first aligned string
        for char in self.firstString:
            string = string + char;
        string = string + "\n";
        # Construct the second aligned string
        for char in self.secondString:
            string = string + char;
        string = string + "\n";
        return string;
    
    # Get the last non-space position of the first string
    def GetLastRealCharOfFirst(self):
        last = 0;
        # Go over all characters and store last non space character
        for index in range(len(self.first)):
            if((self.first)[index] != SPACE):
                last = index;
        return last;
    
    # Get the first non-space position of the first string
    def GetFirstRealCharOfSecond(self):
        first = 0;
        # Go over characters until the first non space character
        for index in range(len(self.first)):
            if((self.second)[index] != SPACE):
                first = index;
                break;
        return first;      

    def GetCommonSubstring(self):
        common = ""
        strike = 0
        strikestart = 0
        maxstrike = 0
        maxstrikestart = 0
        if (self.score != BAD_ALLIGNMENT_SCORE):  
            for index in range(len(self.firstString)):
                firstChar = self.firstString[index]
                secondChar = self.secondString[index]
                if (not firstChar.startswith(MISMATCH_PREFIX) and not firstChar == SPACE and not secondChar == SPACE and not firstChar == JOKER):
                    if (strike == 0):
                        strikestart = index
                    strike += 1
                else:
                    if (strike > maxstrike):
                        maxstrike = strike
                        maxstrikestart = strikestart
                    strike = 0
            if (strike > maxstrike):
                maxstrike = strike
                maxstrikestart = strikestart            
            result = self.firstString[maxstrikestart : maxstrikestart + maxstrike]
            for base in result:
                common += base
            # Try to expand left (try to fix alilgment problem: BLABLABB versus BLABLAB- can also be BLABLABB versus BLABLA-B).
            # This was fixed by changing the order of evaluation of the reoccurence equation, but it can still happen for left extensions
            if (maxstrikestart > 1):
                if (self.firstString[maxstrikestart - 1] == SPACE and self.firstString[maxstrikestart - 2] == self.secondString[maxstrikestart - 1]):
                    common = self.secondString[maxstrikestart - 1] + common
                elif (self.secondString[maxstrikestart - 1] == SPACE and self.secondString[maxstrikestart - 2] == self.firstString[maxstrikestart - 1]):
                    common = self.firstString[maxstrikestart - 1] + common
        return common
        
