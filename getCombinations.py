# -*- getCombinations -*-#
"""
SYNOPSIS
    Function that obtains all possible combinations for the given set of cluster sizes
    
DESCRIPTION
    Function that obtains all possible combinations for the given set of cluster sizes
    
EXAMPLES
    LsCombinations = Get_Combinations(cluster sizes) : the cluster size is formated as [sz0, sz1, ...],
                the output is a library formatted as {sz0: lsCombos0, sz1: lsCombos1},
                if sz_k = sz_l, then the result will be formatted as sz_k = lsCombosk,lsCombosl

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/03
    
"""
from itertools import combinations

def Iter_Get_Combinations(lsCombos, subLsCombos, obs, R, idxLs):
    """
    recursive function to obtain all possible combinations of a cluster or multiple clusters with the same size R 
    input:
        lsCombos - the combination list that contains all possible combinations, structured as [[subLsCombos], ...]
        subLsCombos - the sub-combination list that contains the current combinations, structured as [[combos], ...]
        obs - the list of observations (in int), structured as [0, 1, ...]
        R - size of the cluster
        idxLs - current index of the cluster number, the idx is counted down from the number clusters (w/ same size) to 1
    """
    # create list of combinations from obs into cluster of size R
    currLs0 = [list(combo) for combo in combinations(obs, R)]
    if idxLs == 1:  # if that's the last cluster
        for combo0 in currLs0:
            subLsCombos0 = subLsCombos[:]
            subLsCombos0.append(combo0)
            subLsCombos0.sort(key = lambda c:c[0])  # sort the cluster based on their first observation, *there is no duplicated observation
            if subLsCombos0 not in lsCombos:
                lsCombos.append(subLsCombos0[:])
    else:
        for combo0 in currLs0:
            subLsCombosNext = subLsCombos[:]
            subLsCombosNext.append(combo0)
            obsNext = [e for e in obs if e not in combo0]  # remove the observations in the current cluster from the bank
            Iter_Get_Combinations(lsCombos, subLsCombosNext, obsNext, R, idxLs-1)

def Get_Combinations(clSzs):
    """
    function that obtains all possible combinations for the given set of cluster sizes
    input:
        clSzs - the size of the clusters, structured as [sz0, sz1...]
    return:
        LsCombos - the composite list of combinations, structured as {sz0: lsCombos0, sz1: lsCombos1, ...}
    """
    N_obs = sum(clSzs)
    # compute the frequency where each cluster size occurs
    freqClSz = {}
    for sz in clSzs:
        if sz in freqClSz:
            freqClSz[sz] += 1
        else:
            freqClSz[sz] = 1
    remngN_obs = N_obs
    LsCombos = {}
    currN_combo = 1
    # compute the combination list for each cluster size
    for sz in list(freqClSz):
        obs = range(remngN_obs)  # obervation list is created as [0, 1, ...]
        currLsCombos = []
        currSubLsCombos = []
        Iter_Get_Combinations(currLsCombos, currSubLsCombos, obs, sz, freqClSz[sz])
        LsCombos[sz] = currLsCombos  # add the list that contains all possible combinations for cluster size sz
        currN_combo *= len(currLsCombos)
        remngN_obs -= sz*freqClSz[sz]  # resize and reindex the remaining observations
    return LsCombos
