# -*- clusterTargets -*-#
"""
SYNOPSIS
    Function that reads the newest txt file that contains the SRS target information and compute the optimal clusters
    
DESCRIPTION
    Function that reads the newest txt file that contains the SRS target information and compute the optimal clusters from 1 cluster to up to 3 clusters.
    Both the min sum of square and the min max-dist methods are used.
    
EXAMPLES
    Place the txt file in the same folder of the exe file or one directory up, and the program shall run.

VERSION 0.1
    Add terminal color
AUTHOR
    Becket Hui 2022/04
    
"""
import glob, os, re, sys
import numpy as np
from classTarget import SRSTarget
from classTargetCluster import SRSTargetCluster
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
            print(colored('Computing combination of targets with ', 'yellow') + colored(str(N_cl), 'green') + colored(' cluster.', 'yellow'))
            tgtCl.Calc_Cluster_Sizes(N_cl)
            tgtCl.Calc_Optimal_Cluster()
            print('The optimal cluster is: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for the cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
        else:
            print(colored('Computing optimal ', 'yellow') + colored('Sum of Square', 'green') + colored(' combinations for targets with ', 'yellow') + colored(str(N_cl), 'green') + colored(' cluters.', 'yellow'))
            tgtCl.Calc_Cluster_Sizes(N_cl)
            tgtCl.Calc_Optimal_Cluster()
            print('The optimal clusters are: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for each cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within each of the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
            print(colored('Computing optimal ', 'yellow') + colored('Maximum-Distance to Iso', 'green') + colored(' combinations for targets with ', 'yellow') + colored(str(N_cl), 'green') + colored(' cluters.', 'yellow'))
            tgtCl.Calc_Min_Max_Distance_Cluster()
            print('The optimal clusters are: ')
            print(tgtCl.optimal_cluster)
            print('The ' + colored('iso-center', 'green') + ' for each cluster is')
            print([iso.tolist() for iso in tgtCl.Get_IsoCenters()])
            print('The ' + colored('maximum distance', 'green') + ' from iso within each of the cluster is')
            print(tgtCl.Get_Max_Distance_From_Iso())
        print('')
    input('Hit Enter to Exit...')
    sys.exit()
