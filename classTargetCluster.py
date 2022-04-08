# -*- classTargetCluster -*-#
"""
SYNOPSIS
    A class that includes methods of brute force clustering of SRS targets
    
DESCRIPTION
    The class contains information of the SRS targets including their centroid locations and volume,
    it also contains methods that compute all the combinations of the clusters and find the cluster combination
    with the least sum of square distance, or with the minimum max-cluster to iso distance.
    
EXAMPLES
    srsCluster = SRSTargetCluster(tgts) : initalization of the class
    srsCluster.Calc_Cluster_Sizes(N_cl) : initalize with the number of cluster N_cl
    srsCluster.Calc_Optimal_Cluster() : compute the cluster with the least sum of square distance between targets and isocenter
    srsCluster.Calc_Min_Max_Distance_Cluster() : compute the cluster with the minimum max-target to isocenter distance
    srsCluster.Get_IsoCenters() : get the isocenters of the optimal clusters
    srsCluster.Get_Max_Distance_From_Iso() : get the maximum distance between target and iso within each optimal clusters

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/04
    
"""
import numpy as np
from scipy.spatial import KDTree
from classTarget import SRSTarget
from classNumberOfClusterCombinations import N_Cluster_Combinations
from calcCostFunc import calcCostFunc, calcIsoCenter, calcSqDistFromIso
from getCombinations import Get_Combinations

class SRSTargetCluster:
    """
    SRS Target Cluster class that contains information of the cohort of SRS targets, and fucntion that computes the optimal cluster
    """
    def __init__(self, tgts):
        self.targets = tgts
        self.N_target = len(tgts)
        self.N_cluster = 1
        self.N_cluster_combo = None
        self.max_combo_2calc = 1E5  # keep calculation under 100k
        self.optimal_cluster = []
        self.optimal_cluster_index = []
        self.N_sub_cluster_combo = None
        self.sub_cluster = []
        self.sub_cluster_min_cost = 0
        self.optimal_sub_cluster_index = []

    def Calc_Cluster_Sizes(self, N_cl):
        """
        Calc_Cluster_Sizes computes all the possible cluster sizes and their numbers of combinations
        input:
            SRSTargetCluster object
            N_cl - number of clusters
        return (global):
            N_cluser - number of clusters
            N_cluster_combo - N_Cluster_Combinations object of N_target and N_cl
        """
        if N_cl > self.N_target:
            print('There cannot be more clusters than targets.')
            return
        self.N_cluster = N_cl
        self.N_cluster_combo = N_Cluster_Combinations(self.N_target, self.N_cluster)
        print('There are ' + str(self.N_cluster_combo.N_tot_combinations) + ' possible combinations to form ' + str(self.N_cluster) + ' clusters.')
    
    def Calc_Optimal_Cluster(self):
        """
        Calc_Optimal_Cluster computes the optimal cluster that yields the min cost function
        input:
            SRSTargetCluster object
        """
        if self.N_cluster_combo is None:
            print('Please run Calc_Cluster_Sizes first before finding optimal cluster.')
            return
        # first create the np array of the sub-cluster that contains [idx, ctr x, ctr y, ctr z, vol fctr]
        subClNpArr = self.Pre_Cluster_Targets()
        self.sub_cluster_min_cost = calcCostFunc(subClNpArr[:, 1:4], subClNpArr[:, 4]) + 1  # maximum cost function is single cluster
        for clSzs in self.N_sub_cluster_combo.Ls_cluster_size:
            lsCombos = Get_Combinations(clSzs.sizes)
            self.Iter_Calc_Cost_Function(subClNpArr, lsCombos, 0, 0.0, [])
            del lsCombos  # free up space
        # uncluster the sub-cluster and convert cluster of sub-cluster to cluster of targets
        self.Unwrap_Sub_Cluster()

    def Iter_Calc_Cost_Function(self, subClNpArr, lsCombos, idxClSz, currCost_0, currCluster_0):
        """
        Iter_Calc_Cost_Function is a recursive funtion that computes the min cost function from the input cluster combinations
        input:
            SRSTargetCluster object
            subClNpArr - sub-cluster numpy array, structured as [[idx, ctr x, ctr y, ctr z, vol fctr], ...]
            lsCombos - list of cluster combinations, structured as [[idx1, idx2, ...], ...]
            idxClSz - index of cluster size, go from 0 to len(lsCombos)-1
            currCost_0 - current cost funtion composited from all the previous cluster size indices
            currCluster_0 - current formed cluster from the previous cluster size indices, structured as [[idxA, idxB, ...], ...]
        return (global):
            sub_cluster_min_cost - will be updated if new minimum is found
            optimal_sub_cluster_index - will be the new currCluster if a new minimum is found
        """
        for combo in lsCombos[list(lsCombos)[idxClSz]]:  # loop over all combinations for the cluster size list(lsCombos)[idxClSz] as lsCombos = [combo, combo, ...]
            currCost = currCost_0
            currCluster = currCluster_0[:]
            for subcombo in combo:  # subcombo(s) within a combo have same size, loop over all subcombinations as combo = [[subcombo], [subcombo], ...]
                currCost = currCost + calcCostFunc(subClNpArr[subcombo, 1:4], subClNpArr[subcombo, 4])
                currCluster.append(subClNpArr[subcombo, 0])
            if currCost < self.sub_cluster_min_cost:  # if the current cost is smaller than the minimum cost, proceed, otherwise go to next loop
                if idxClSz == len(lsCombos)-1:  # if that's the last cluster size
                    self.sub_cluster_min_cost = currCost
                    self.optimal_sub_cluster_index = currCluster
                    # print('The cost function is ' + str(currCost) + ' for cluster:')  # DEBUG
                    # print([clNp.tolist() for clNp in currCluster])  # DEBUG
                else:
                    subClNpArr_Next = np.delete(subClNpArr, (sum(combo, [])), axis = 0)  #remove the combinations already computed and move to next loop
                    self.Iter_Calc_Cost_Function(subClNpArr_Next, lsCombos, idxClSz+1, currCost, currCluster)

    def Calc_Min_Max_Distance_Cluster(self):
        """
        Min_Max_Distance_Cluster computes the cluster that yields the min(max distance from isocenter)
        input:
            SRSTargetCluster object
        """
        if self.N_cluster_combo is None:
            print('Please run Calc_Cluster_Sizes first before finding optimal cluster.')
            return
        # first create the np array of the sub-cluster that contains [idx, ctr x, ctr y, ctr z, vol fctr]
        subClNpArr = self.Pre_Cluster_Targets()
        # create the min cost array as the [(max dist)^2 + 1] x N_cluster
        self.sub_cluster_min_cost = [max(calcSqDistFromIso(subClNpArr[:, 1:4])) + 1] * self.N_cluster
        for clSzs in self.N_sub_cluster_combo.Ls_cluster_size:
            lsCombos = Get_Combinations(clSzs.sizes)
            self.Iter_Calc_Max_Dist(subClNpArr, lsCombos, 0, [], [])
            del lsCombos  # free up space
        # uncluster the sub-cluster and convert cluster of sub-cluster to cluster of targets
        self.Unwrap_Sub_Cluster()

    def Iter_Calc_Max_Dist(self, subClNpArr, lsCombos, idxClSz, currMaxDist2_0, currCluster_0):
        """
        Iter_Calc_Cost_Function is a recursive funtion that computes the min cost function from the input cluster combinations
        input:
            SRSTargetCluster object
            subClNpArr - sub-cluster numpy array, structured as [[idx, ctr x, ctr y, ctr z, vol fctr], ...]
            lsCombos - list of cluster combinations, structured as [[idx1, idx2, ...], ...]
            idxClSz - index of cluster size, go from 0 to len(lsCombos)-1
            currMaxDist2_0 - current max distance from iso for all the previous cluster size indices
            currCluster_0 - current formed cluster from the previous cluster size indices, structured as [[idxA, idxB, ...], ...]
        return (global):
            sub_cluster_min_cost - will be updated if new minimum is found
            optimal_sub_cluster_index - will be the new currCluster if a new minimum is found
        """
        for combo in lsCombos[list(lsCombos)[idxClSz]]:  # loop over all combinations for the cluster size list(lsCombos)[idxClSz] as lsCombos = [combo, combo, ...]
            currMaxDist2 = currMaxDist2_0[:]
            currCluster = currCluster_0[:]
            for subcombo in combo:  # subcombo(s) within a combo have same size, loop over all subcombinations as combo = [[subcombo], [subcombo], ...]
                currMaxDist2.append(max(calcSqDistFromIso(subClNpArr[subcombo, 1:4])))
                currCluster.append(subClNpArr[subcombo, 0])
            # If the current max distance within clusters is smaller than the min(max distance), proceed, otherwise go to next loop
            if self.Max_Dist_1_LE_Max_Dist_2(currMaxDist2, self.sub_cluster_min_cost):
                if idxClSz == len(lsCombos)-1:  # if that's the last cluster size
                    self.sub_cluster_min_cost = currMaxDist2[:]  # the comparison function has sorted the distance list
                    self.optimal_sub_cluster_index = currCluster
                    # print('The max distances from iso are ' + str(np.sqrt(currMaxDist2)) + ' for cluster:')  # DEBUG
                    # print(currCluster)  # DEBUG
                else:
                    subClNpArr_Next = np.delete(subClNpArr, (sum(combo, [])), axis = 0)  #remove the combinations already computed and move to next loop
                    self.Iter_Calc_Max_Dist(subClNpArr_Next, lsCombos, idxClSz+1, currMaxDist2, currCluster)
       
    def Pre_Cluster_Targets(self):
        """
        Pre_Cluster_Targets pre cluster the targets into sub-clusters before calculating the optimal cluster
        input:
            SRSTargetCluster object
        return (global):
            sub_cluster - the list of subclusters in target indices
        """
        if self.N_cluster_combo is None:
            print('Please run Calc_Cluster_Sizes first before running pre-clustering.')
            return
        # first compute the number of sub-clusters that target set needs to be reduced to
        N_subCl = self.Get_N_Targets_4_Calc()
        if N_subCl < self.N_target:
            print('The targets will be pre-clustered to ' + str(N_subCl) + ' clusters before computing optimal clusters.')
        # start pre-clustering targets
        subClLs = [[idxTgt] for idxTgt in range(self.N_target)]
        subClCtrArr = np.empty((self.N_target, 3))
        subClVolArr = np.empty((self.N_target))
        for idxTgt in range(self.N_target):
            subClCtrArr[idxTgt,:] = self.targets[idxTgt].centroid
            subClVolArr[idxTgt] = self.targets[idxTgt].volume
        while len(subClLs) > N_subCl:
            # use K-D Tree to solve for nearest neighbour
            # pre-group the nearest neighbours one-by-one until number of sub-clusters equals N_subCl
            tree = KDTree(subClCtrArr)
            dists, idcs = tree.query(subClCtrArr, k=2)
            idxMin = np.argmin(dists[:,1])
            # culster the two points with the smallest distance
            subClLs, subClCtrArr, subClVolArr = self.Merge_Two_Points(subClLs, subClCtrArr, subClVolArr, idcs[idxMin,:])
            # print('The sub-cluster indices are:')  #DEBUG
            # print(subClLs)  #DEBUG
            # print('The average centroid positions of them are:')  #DEBUG
            # print(subClCtrArr)  #DEBUG
        subClVolFctrArr = np.array([len(sc) for sc in subClLs])  # use the number of targets in subcluster as volume factor
        # assign precluster and create np array for optimal cluster calculation
        self.sub_cluster = subClLs
        subClNpArr = np.empty((N_subCl, 5))
        for idxSubCl in range(N_subCl):
            subClNpArr[idxSubCl, :] = np.concatenate(([idxSubCl], subClCtrArr[idxSubCl, :], [subClVolFctrArr[idxSubCl]]))
        return subClNpArr

    def Get_N_Targets_4_Calc(self):
        """
        Compute the max number of targets that Calc_Optimal_Cluster will handle, such that the resulting number of combinations will be less than max_combo_2calc
        input:
            SRSTargetCluster object
        return:
            currN_tgt - max number of allowed sub-clusters/points for optimal cluster calculation
            N_sub_cluster_combo (global) - N_Cluster_Combinations object with N_subCl instead of N_target
        """
        currN_tgt = self.N_target
        nextN_tgt = max(currN_tgt//2, self.N_cluster)
        N_combos = self.N_cluster_combo.N_tot_combinations
        while N_combos > self.max_combo_2calc:  # if number of combinations is more than the max to compute
            N_clCombo = N_Cluster_Combinations(nextN_tgt, self.N_cluster)
            N_combos = N_clCombo.N_tot_combinations
            currN_tgt = nextN_tgt
            nextN_tgt = max(currN_tgt//2, self.N_cluster)
            # print('For number of targets = ' +  str(currN_tgt) + ', number of combinations = ' + str(N_combos) + '.')  # DEBUG
            if currN_tgt == nextN_tgt:  #repeated loop, break
                break
        # after while loop, the number of combinations is less thatn max to compute for currN_tgt
        # next increase currN_tgt to find the max currN_tgt that produces the number of combinations less than max to compute
        if currN_tgt < (self.N_target - 1):  # only run if currN_tgt is not N_target or N_target-1
            nextN_tgt = currN_tgt
            while N_combos < self.max_combo_2calc:  # if number of combinations for nextN_tgt is over max to compute, then curr_Ntgt is our currN_tgt is the max 
                nextN_tgt += 1
                currN_tgt = nextN_tgt - 1
                if nextN_tgt == self.N_target:  # if the next tested target is the total number of target
                    break
                N_clCombo = N_Cluster_Combinations(nextN_tgt, self.N_cluster)
                N_combos = N_clCombo.N_tot_combinations
                # print('For number of targets = ' +  str(nextN_tgt) + ', number of combinations = ' + str(N_combos) + '.')  # DEBUG
        self.N_sub_cluster_combo = N_Cluster_Combinations(currN_tgt, self.N_cluster)
        return currN_tgt  # number of maximum allowed targets

    def Merge_Two_Points(self, subClLs, ctrArr, volArr, idxMerge):
        """
        Merge_Two_Points function merge two points/sub-clusters to form a new redcued dimension sub-clusters, it recalculates the centroid and volume of the new sub-cluster
        input:
            SRSTargetCluster object
            subClLs - the sub-cluster list that contains the target indices of the sub-clusters 
            ctrArr - the centroid array of the sub-clusters
            volArr - the volume array of the sub-clusters
            idxMerge - the list of two indices of the sub-cluster list to merge
        return:
            subClLs - the new sub-cluster list that contains the target indices of the sub-clusters 
            ctrArr - the centroid array of the new sub-clusters
            volArr - the volume array of the new sub-clusters 
        """
        # first put points from idxMerge[1] to the sub-cluster list subClLs[idxMerge[0]]
        for tgt in subClLs[idxMerge[1]]:
            subClLs[idxMerge[0]].append(tgt)
        del subClLs[idxMerge[1]]
        # then recompute the centroid location of the newly merged sub-cluster
        subClCtrArr = np.empty((len(subClLs[idxMerge[0]]),3))
        for idxTgt in range(len(subClLs[idxMerge[0]])):
            subClCtrArr[idxTgt,:] = self.targets[subClLs[idxMerge[0]][idxTgt]].centroid
        ctrArr[idxMerge[0],:] = np.mean(subClCtrArr, axis=0)  # average centroid of targets in sub-cluster
        ctrArr = np.delete(ctrArr, idxMerge[1], 0)
        # finally compute the effective volume of the newly merged cluster
        invSum = 0
        for idxTgt in range(len(subClLs[idxMerge[0]])):
            invSum = 1/(self.targets[subClLs[idxMerge[0]][idxTgt]].volume)
        volArr[idxMerge[0]] = 1/invSum
        volArr = np.delete(volArr, idxMerge[1], 0)
        return subClLs, ctrArr, volArr

    def Unwrap_Sub_Cluster(self):
        """
        Unwrap_Sub_Cluster unclusters the subcluster and expands the clusters of sub-clusters back to clusters of targets
        input:
            SRSTargetCluster object
        return (global):
            optimal_cluster_index - the optimal clusters in target indices
            optimal_cluster - the optimal clusters in target names
        """
        if self.sub_cluster_min_cost == 0:
            print('Please run Calc_Optimal_Cluster before unwrapping sub clusters.')
            return
        self.optimal_sub_cluster_index = [cl.astype(int).tolist() for cl in self.optimal_sub_cluster_index]  # cast the index from float to int
        opClIdx = []
        opClName = []
        for subCl in self.optimal_sub_cluster_index:
            currClIdx = []
            for clIdx in subCl:
                currClIdx.extend(self.sub_cluster[clIdx])
            currClIdx.sort()
            opClIdx.append(currClIdx)
            opClName.append([self.targets[idx].name for idx in currClIdx])
        self.optimal_cluster_index = opClIdx
        self.optimal_cluster = opClName

    def Max_Dist_1_LE_Max_Dist_2(self, maxDistLs1, maxDistLs2):
        """
        Max_Dist_1_LE_Max_Dist_2 compares Ls1 with Ls2 and return True if the largest element in Ls1 is less than those in Ls2,
        it also return True if all elements in Ls1 are the same as the top elements in Ls2.
        It is used only within Iter_Calc_Max_Dist.
        """
        maxDistLs1.sort(reverse = True)
        maxDistLs2.sort(reverse = True)
        # comparison will only compare the maximum distance from the lists
        for idxCl in range(len(maxDistLs1)):
            if maxDistLs1[idxCl] < maxDistLs2[idxCl]:
                return True
            if maxDistLs1[idxCl] > maxDistLs2[idxCl]:
                return False
            # if the distances in the list are equal, the comparison will go to the next largest distances
        return True  # if all distances are equal, return true

    def Get_IsoCenters(self):
        """
        Compute the group of isocenters for the optimal clusters
        """
        isoCtrs = []
        if not self.optimal_cluster_index:
            print('Please run Clac_Optimal_Cluster first before finding the cluster isocenter.')
        else:
            for cl in self.optimal_cluster_index:
                ctrArr = np.empty((len(cl),3))
                for idxTgt in range(len(cl)):
                    ctrArr[idxTgt,:] = self.targets[cl[idxTgt]].centroid
                isoCtrs.append(calcIsoCenter(ctrArr))
        return isoCtrs

    def Get_Max_Distance_From_Iso(self):
        """
        Compute the max distance between the target and the isocenter within each cluster
        """
        maxDist = []
        if not self.optimal_cluster_index:
            print('Please run Clac_Optimal_Cluster first before finding the distance from isocenter.')
        else:
            for cl in self.optimal_cluster_index:
                ctrArr = np.empty((len(cl),3))
                for idxTgt in range((len(cl))):
                    ctrArr[idxTgt, :] = self.targets[cl[idxTgt]].centroid
                sqDist = calcSqDistFromIso(ctrArr)
                maxDist.append(np.sqrt(max(sqDist)))
        return maxDist

    # def Get_Current_Cost_Function(self):  #DEBUG
    #     costFunc = 0.0
    #     if not self.optimal_cluster_index:
    #         print('Please run Clac_Optimal_Cluster first before finding the cost function.')
    #     else:
    #         for cl in self.optimal_cluster_index:
    #             ctrArr = np.empty((len(cl),3))
    #             volFctrArr = np.empty((len(cl)))
    #             for idxTgt in range(len(cl)):
    #                 ctrArr[idxTgt, :] = self.targets[cl[idxTgt]].centroid
    #                 volFctrArr[idxTgt] = self.targets[cl[idxTgt]].volume_factor
    #             costFunc += calcCostFunc(ctrArr, volFctrArr)
    #     return costFunc