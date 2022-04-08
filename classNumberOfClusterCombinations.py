# -*- classNumberOfClusterCombinations -*-#
"""
SYNOPSIS
    Class that compute the number of combinations given the number of targets and clusters
    
DESCRIPTION
    Class that contains methods to compute the number of combinations given the number of targets and clusters
    
EXAMPLES
    N_Cluster_Combo = N_Cluster_Combinations(N_targets, N_cluster)

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/03
    
"""
from calcPermutationAndCombination import calcCombinations as nCr, calcPermutaions as nPr
from math import factorial

class N_Cluster_Combinations:
    """
    N_Cluster_Combinations class to get number of possible combinations for N_obs observations partitioned into N_cl clusters
    """
    def __init__(self, N_obs, N_cl):
        # check input types
        if not self.Chk_Int(N_obs):
            print('Input number of observations is not positive integer.\n')
            return
        if not self.Chk_Int(N_cl):
            print('Input number of clusters is not positive integer.\n')
            return
        if N_cl > N_obs:
            print('Input number of clusters is more than number of observations.\n')
            return
        self.N_observation = N_obs
        self.N_cluster = N_cl
        # create all possible cluster sizes and calculate the corresponding number of combinations
        self.Ls_cluster_size = []
        LsClSz = self.Create_Cluster_Sizes()
        for idxCl in range(len(LsClSz)):
            self.Ls_cluster_size.append(self.Cluster_Size(self, LsClSz[idxCl]))
        self.N_tot_combinations = sum(cl.N_combo for cl in self.Ls_cluster_size)

    def Create_Cluster_Sizes(self):
        """
        function that creates the all combinations of dimensions in each cluster
        input:
            N_Cluster_Combinations - the N_Cluster_Combinations class object
        return:
            LsClSz - list of dimension combinations
        """
        # initiate range of dimension for the first cluster
        minClSz_0 = self.N_observation // self.N_cluster
        maxClSz_0 = self.N_observation - self.N_cluster + 1
        # compute all possible number of members in each cluster
        currClSzs = [0]*self.N_cluster
        LsClSz = []
        self.Iter_Calc_N_Members(LsClSz, currClSzs, self.N_cluster, minClSz_0, maxClSz_0)
        return LsClSz
  
    def Iter_Calc_N_Members(self, LsClSz, currClSzs, N_remngCl, minClSz, maxClSz):
        """
        recursive function to calculate the dimension combinations in each cluster
        input:
            N_Cluster_Combinations - the N_Cluster_Combinations class object
            LsClSz - list that contains all possible dimension combinations in each cluster, structured as [currClSz, ...]
            currClSz - current dimension in each cluster, structured as [dim0, ...]
            N_remngCl - index of the cluster, counting from number of clusters to 1
            minClSz - minimum dimension of the current cluster
            maxClSz - maximum dimension of the current cluster
        """
        if N_remngCl == 1:  # if this is the last cluster
            currClSzs[self.N_cluster - 1] = self.N_observation - sum(currClSzs)
            currClSzs.sort(reverse = True)  # sort the cluster sizes from large to small
            if currClSzs not in LsClSz and min(currClSzs) > 0:
                LsClSz.append(currClSzs[:])  # creates a copy for currClSzs and append to Ls
        else:  # if there are other remaining cluster to calculate number of combination
            for dim in range(maxClSz, minClSz-1, -1):
                currClSzs[self.N_cluster-N_remngCl] = dim
                # reset dimension of subsequent clusters
                for idxDim in range(self.N_cluster-N_remngCl + 1, self.N_cluster): currClSzs[idxDim] = 0
                N_remngCl_Next = N_remngCl-1
                minClSz_Next = (self.N_observation - sum(currClSzs)) // N_remngCl_Next
                maxClSz_Next = min(self.N_observation - sum(currClSzs) - N_remngCl_Next + 1, dim)
                # compute dimension for the next cluster
                if maxClSz_Next >= minClSz_Next:
                    self.Iter_Calc_N_Members(LsClSz, currClSzs, N_remngCl_Next, minClSz_Next, maxClSz_Next)

    class Cluster_Size:
        """
        Cluster_Size class that contains the cluster sizes and its number of combinations
        """
        def __init__(self, Obj_NClCombo, clSzs):
            self.sizes = clSzs
            if Obj_NClCombo.N_observation == sum(self.sizes):
                self.N_combo = self.Calc_N_Combo(Obj_NClCombo.N_observation)
            else:
                print('Combined clusters dimension is not equal to number of observations.\n')
                self.N_combo = 0

        def Calc_N_Combo(self, N_Obs):
            """
            function that computes the number of combinations in the given cluster dimension
            input:
                Cluster_Size - the Cluster_Size class object
                N_Obs - number of observations
            return:
                N_combo - number of combinations
            """
            N_cluster = len(self.sizes)
            N_remngObs = N_Obs
            N_combo = 1
            # compute combinations for each cluster and multiply
            for idxCl in range(N_cluster-1):
                N_clMem = self.sizes[idxCl]
                N_combo *= nCr(N_remngObs, N_clMem)
                N_remngObs -= N_clMem
            # remove permutation between clusters with same dimension
            freqClSz = {}
            for dimCl in self.sizes:
                if dimCl in freqClSz:
                    freqClSz[dimCl] += 1
                else:
                    freqClSz[dimCl] = 1
            N_clPermute = 1
            for freq in freqClSz.values():
                N_clPermute *= factorial(freq)
            return int(N_combo/N_clPermute)

    def Chk_Int(self, x):
        # check if input x is integer or list of integers
        xIsInt = None
        try:
            xIsInt = False
            if isinstance(x, list):
                if all(isinstance(xi, int) for xi in x):
                    if all((xi > 0) for xi in x): xIsInt = True
            if isinstance(x, int):
                if x > 0: xIsInt = True
            return xIsInt
        except:
            print('Cannot determine if input is integer.\n')
            return xIsInt


