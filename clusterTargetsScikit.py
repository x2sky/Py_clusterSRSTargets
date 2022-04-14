# -*- clusterTargetsScikit -*-#
"""
SYNOPSIS
    Function that reads the newest txt file that contains the SRS target information and compute the optimal clusters
    
DESCRIPTION
    Function that reads the newest txt file that contains the SRS target information and compute the optimal clusters from 1 cluster to up to 3 clusters.
    The l1, l2 hierarchial and k-means clustering methods are used
    
EXAMPLES
    Place the txt file in the same folder of the exe file or one directory up, and the program shall run.

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/04
    
"""
import glob, os, re, sys
import numpy as np
from classTarget import SRSTarget
from classScikitCluster import SRSTargetCluster
from termcolor import colored

if __name__ == '__main__':
    # first read file and extract the targets' info
    currDir = os.getcwd()
    # currDir = r'M:\Eclipse\Scripting\SRSCluster'
    fileLs = glob.glob(os.path.join(currDir, '*_srstgt.txt'))
    if not fileLs:
        currDir = os.path.dirname(currDir)  # move one directory up
        fileLs = glob.glob(os.path.join(currDir, '*_srstgt.txt'))
    if not fileLs:  # if still cannot find file, exit
        input('Cannot find any srstgt.txt file, hit Enter to Exit...')
        sys.exit()
    # read the latest file
    lastestFile = max(fileLs, key = os.path.getmtime)
    filePtr = open(lastestFile, 'r')
    srsTargets = []
    for txt in filePtr:
        info = txt.split(', ')
        srsTargets.append(SRSTarget(info[0], np.array(info[1:4]).astype(float), float(info[4])))
    filePtr.close()
    # perform clustering
    os.system('mode con: cols=250')
    m = re.search('(1100[0-9]+)_', lastestFile)
    if m:
        print('Start target clustering for patient ' + colored(m.group(1), 'yellow') + '...')
    else:
        print('Start target clustering...')
    tgtCl = SRSTargetCluster(srsTargets)
    for N_cl in range(1,min(4, tgtCl.N_target+1),1):  # calculate from 1 up to 3 or N_target clusters
        if N_cl == 1:
            print(colored('Computing combination of targets with ' + str(N_cl) + ' cluster.', 'yellow'))
            tgtCl.Calc_Optimal_Cluster(N_cl, 'kmeans')
            print('The optimal cluster is: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for the cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
            print('')
        else:
            print(colored('Computing optimal combinations for targets with ' + str(N_cl) + ' cluters.', 'yellow'))
            tgtCl.Calc_Optimal_Cluster(N_cl, 'hrchy2')
            print('Using ' + colored('l2 hierarchical', 'green') + ' clustering, the optimal clusters are: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for each cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within each of the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
            print('')

            tgtCl.Calc_Optimal_Cluster(N_cl, 'hrchy1')
            print('Using ' + colored('l1 hierarchical', 'green') + ' clustering, the optimal clusters are: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for each cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within each of the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
            print('')

            tgtCl.Calc_Optimal_Cluster(N_cl, 'kmeans')
            print('Using ' + colored('k-means hierarchical', 'green') + ' clustering, the optimal clusters are: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for each cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within each of the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
            print('')
    input('Hit Enter to Exit...')
    sys.exit()
