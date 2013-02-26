#!/usr/bin/env python
__version__ = '1.0'
__date__ = 'May 18, 2008'
__author__ = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI'
__url__ = 'http://quwubin.wordpress.com.cn/'
__license__ = 'GPL'

import LightSeq as LS

ambiguous_dna_complement = {
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "A": "T", 
    "C": "G", 
    "G": "C", 
    "T": "A", 
    "M": "K", 
    "m": "k",
    "R": "Y", 
    "r": "y",
    "W": "W", 
    "w": "w",
    "S": "S", 
    "s": "s",
    "Y": "R", 
    "y": "r",
    "K": "M", 
    "k": "m",
    "V": "B", 
    "v": "b",
    "H": "D", 
    "h": "d",
    "D": "H", 
    "d": "h",
    "B": "V", 
    "b": "v",
    "X": "X", 
    "x": "x",
    "N": "N", 
    "n": "n",
    "-": "-", 
}

def comSeq(s):
    tmp = ''
    for i in range(len(s)):
        tmp = tmp + ambiguous_dna_complement[s[i]]

    return tmp

def revSeq(s):
    """Reverses a string given to it."""
    return s[::-1]

def reverseString(s):
    """Reverses a string given to it."""
    return s[::-1]

def revComSeq(s):
    s = revSeq(s)
    s = comSeq(s)
    return s

class seq():
    def __init__(self, sequence):
        self.reverse = revSeq(sequence)
        self.complement = comSeq(sequence)
        self.rev_com = revComSeq(sequence)


# The master test function
def test():
    seq = 'CGTTGA'
    print LS.seq(seq).complement
    print '*' * 80
    print seq
    
if __name__ == "__main__":
    test()
    print "Done."
