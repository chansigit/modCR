from time import *

HASH_EXP = 17
MOD = 2**HASH_EXP - 1
#BASE = 4
BASE = 5
BASE = 6
#BASE = 7
#ALPHABET_VALUES = {"A":0, "G":1 , "C":2, "T":3} sign(AACT) = sign(CT)
ALPHABET_VALUES = {"A":1, "G":2 , "C":3, "T":4}
ALPHABET_VALUES = {"A":1, "G":2 , "C":3, "T":4, "N":5}
#ALPHABET_VALUES = {"A":1, "G":2 , "C":3, "T":4, "N":5, "D":6}
PREFIX = 0
SUFFIX = 1

# General method for calculating hash values of prefixes/suffixes
def CalcHashValues(strindex, string, minoverlap, side, signatures, m, b = BASE):
    if (side == PREFIX):
        return CalcPrefixValues(strindex, string, minoverlap, signatures, m, b)
    elif (side == SUFFIX):
        return CalcSuffixValues(strindex, string, minoverlap, signatures, m, b)
    else:
        assert(false)
        
# Methods using regular (non incremental) hash function
def CalcHPrefixValues(strindex, string, minoverlap):
    results = []
    for index in range(minoverlap, len(string), 1):
        results.append(hash(string[:index]))
    return results

def CalcHSuffixValues(strindex, string, minoverlap):
    results = []
    for index in range(minoverlap, len(string), 1):
        results.append(hash(string[-index:]))
    return results

def CalcKRValue(string, m, b=BASE):
    result = 0
    l = len(string)
    for index in range(len(string)):
        result = (result + (ALPHABET_VALUES[string[index]] * (b ** (l-index-1)))) % m
    return result
    
# Methods using incremental hash function
def CalcPrefixValues(strindex, string, minoverlap, signatures, m, b):
    results = []
    temp = 0
    minimalsize = minoverlap - 1
    length = len(string)
    for index in range(length):
        current = ALPHABET_VALUES[string[index]]
        temp = ((temp * b) + current) % m
        if (index >=  minimalsize):
            if (index < length - 1):
                results.append(temp)
            else:
                if (signatures != None):
                    signatures[temp] = signatures.get(temp, set())
                    signatures[temp].add(strindex)
    return results

def CalcSuffixValues(strindex, string, minoverlap, signatures, m, b):
    results = []
    temp = 0
    exp = 1
    minimalsize = minoverlap - 1
    for index in range(len(string) - 1):
        current = ALPHABET_VALUES[string[-(index+1)]]
        temp = ((exp * current) + temp) % m
        exp = (exp * b) % m
        if (index >=  minimalsize):
            results.append(temp)
    return results

# Extending incremental hash functions
# An incremental hash is necessary here to have no need in underlying string.
# Return value of the method is a mapping of extending symbols to hash values.
def ExtendKR(value, multiplier, side, m, b = BASE):
    if (side == PREFIX):
        return ExtendKRLeft(value, multiplier, m, b)
    elif (side == SUFFIX):
        return ExtendKRRight(value, multiplier, m, b)
    else:
        assert(false)
        
def ExtendKRRight(value, multiplier, m, b):
    extensions = dict()
    temp = (value * b) % m
    for symbol in ALPHABET_VALUES.keys():
        tempValue = ((temp + ALPHABET_VALUES[symbol]) % m)
        extensions[symbol] = tempValue
    return extensions

def ExtendKRLeft(value, multiplier, m, b):
    extensions = dict() 
    for symbol in ALPHABET_VALUES.keys():
        tempValue = ((multiplier * ALPHABET_VALUES[symbol] + value) % m)
        extensions[symbol] = tempValue    
    return extensions

def VerifyOverlap(left, right, mo):
    for index in range(mo, min(len(left), len(right))):
        if (left[-index:] == right[:index]):
            return True
    return False

def VerifyCertainOverlap(left, right, overlap, extensionLength, overlaper):
    t1 = clock()
    if (overlap + extensionLength != len(right)):
        t2 = clock()
        #overlaper.vtimer += (t2-t1)
        #overlaper.vc += 1
        overlaper.falseVerLength += 1
        return False
    if (left[-overlap:] == right or right[:overlap] == left):
        print ("Containment!")
        t2 = clock()
        overlaper.vtimer += (t2-t1)        
        overlaper.vc += 1
        overlaper.falseVerBases += 1
        return False
    temp = 0
    length = len(left)
    for index in range(overlap):
        temp += 1
        if (left[length-overlap+index] != right[index]):
            t2 = clock()
            overlaper.vtimer += (t2-t1)    
            overlaper.vc += 1
            overlaper.falseVerBases += 1
            overlaper.checkedbases += temp
            return False
    t2 = clock()
    overlaper.vtimer += (t2-t1)
    overlaper.trueVer += 1
    return True    


def IsInnerCollision(strindex, string, minoverlap, signatures, m, b):
    prefixes = CalcPrefixValues(strindex, string, minoverlap, signatures, m, b)
    for val in signatures:
        if val in prefixes:
            return True
    return False

def InterCollisions(reads, minoverlap, m, b):
    signatures = {}
    for index in range(len(reads)):
        CalcPrefixValues(index, reads[index], minoverlap, signatures, m, b)
    s = {}
    for index in range(len(reads)):
        prefixes = CalcPrefixValues(index, reads[index], minoverlap, s, m, b)        
        for prefix in prefixes:
            if prefix in signatures:
                others = signatures[prefix]
                for other in others:
                    overlap = 0
                    for os in range(minoverlap, 51):
                        if (reads[index][:os] == reads[other][:os]):
                            overlap = os
                    if overlap > 0:
                        print ("Inter Collision", index, other, overlap)
                            
                
    
s1 = "AGTTGACGTA"
s2 = "AAGCTGTTTT"
t1 = clock()
for index in range(1):
    res = CalcSuffixValues(0, s1, 15, {}, MOD, BASE)
t2 = clock()
#print (t2-t1)
#print (CalcPrefixValues(0, s1, 2, {}, 2 ** 29 - 1, 6))
#print (CalcPrefixValues(1, s2, 2, {}, 2 ** 29 - 1, 6))
#print (CalcSuffixValues(0, s1, 2, {}, 2 ** 29 - 1, 6))
#print (CalcSuffixValues(1, s2, 2, {}, 2 ** 29 - 1, 6))
