# -*- calcPermutationAndCombination -*-#
"""
SYNOPSIS
    Functions that compute nCr and nPr
    
DESCRIPTION
    Functions that compute nCr and nPr
    
EXAMPLES
    calcCombinations(n, r) : nCr
    calcPermutations(n, r) : nPr

VERSION 0.0
    Initial version
AUTHOR
    Becket Hui 2022/03
    
"""
import operator as op
from functools import reduce

def calcCombinations(n, r):
    """
    Choose r elements from a set of n elements
    """
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom

def calcPermutaions(n, r):
    """
    Pick r element sequence from a set of n elements
    """
    result = 1
    for i in range(n, n-r, -1):
        result *= i
    return result