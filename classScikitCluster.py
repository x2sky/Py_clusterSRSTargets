# -*- classScikitCluster -*-#
"""
SYNOPSIS
    A class that includes methods of scikit learn clustering of SRS targets
    
DESCRIPTION
    The class contains information of the SRS targets including their centroid locations and volume,
    it also contains methods that find the cluster using K-means and hierarchial clustering
    
EXAMPLES
    srsCluster = SRSTargetCluster(tgts) : initalization of the class
    srsCluster.Calc_Optimal_Cluster(N_cl, algorithm) : cluster the targets into N_cl clusters using the given algorithm
    srsCluster.Get_IsoCenters() : get the isocenters of the optimal clusters
    srsCluster.Get_Max_Distance_From_Iso() : get the maximum distance between target and iso within each optimal clusters

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/04
    
"""
import numpy as np
from classTarget import SRSTarget
from calcCostFunc import calcIsoCenter, calcSqDistFromIso
from sklearn.cluster import KMeans, AgglomerativeClustering

class SRSTargetCluster:
    """
    SRS Target Cluster class that contains information of the cohort of SRS targets, and fucntion that computes the optimal cluster
    """
    def __init__(self, tgts):
        self.targets = tgts
        self.N_target = len(tgts)
        self.N_cluster = 1
        self.optimal_cluster = []
        self.optimal_cluster_index = []

    def Calc_Optimal_Cluster(self, N_cl, algorithm = 'hrchy2'):
        """
        """
        if not isinstance(N_cl, int):
            print('Input number of cluster is not an integer, exiting...')
            return
        self.N_cluster = min(self.N_target, max(1, N_cl))
        # setup numpy array for clustering, the format structured as [[target idx, centroid coordinates x3, volume factor], ...]
        tgtNpArr = np.empty((self.N_target,5))
        for idxTgt in range(self.N_target):
            tgtNpArr[idxTgt, :] = np.concatenate(([idxTgt], self.targets[idxTgt].centroid, [self.targets[idxTgt].volume_factor]))
        # start clustering
        if algorithm == 'hrchy2':
            clResult = AgglomerativeClustering(n_clusters=N_cl, linkage='complete', compute_full_tree='auto', affinity='l2').fit(tgtNpArr[:, 1:4])
        elif algorithm == 'hrchy1':
            clResult = AgglomerativeClustering(n_clusters=N_cl, linkage='complete', compute_full_tree='auto', affinity='l1').fit(tgtNpArr[:, 1:4])
        elif algorithm == 'kmeans':
            clResult = KMeans(n_clusters=N_cl, init='k-means++', algorithm='full', n_init=20, max_iter=10000).fit(tgtNpArr[:, 1:4])
        else:
            print('Input algorithm is not recognized, exiting...')
            return
        # assgin cluster result
        self.optimal_cluster = []
        self.optimal_cluster_index = []
        for idxCl in range(self.N_cluster):
            cluster0 = (np.where(clResult.labels_ == idxCl)[0]).tolist()  # the clustered indices in tgtNpArr
            cluster = (tgtNpArr[cluster0, 0].astype(int)).tolist()  # converted to indices in self.targets
            self.optimal_cluster_index.append(cluster)
            self.optimal_cluster.append([self.targets[idx].name for idx in cluster])
        return

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