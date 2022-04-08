# -*- calcCostFunc -*-#
"""
SYNOPSIS
    Functions that compute the cost function, square distance, and iso-center
    
DESCRIPTION
    These functions are used by classTargetCluster to compute the sum of square distance cost function;
    the square distance from isocenter for a set of target postions; and the optimal isocenter
    
EXAMPLES
    calcCostFunc(centroidArray, volFactorArray) : compute the sum of square cost function
    calcSqDistFromIso(centroidArray, isoCenter (optional)) : compute the square distances of targets from iso,
                    if no iso provided, will use calcIsoCenter to compute the optimal iso
    calcIsoCenter(centroidArray) : compute the isocenter using the min diameter sphere method

VERSION 0.0
    First functional version
AUTHOR
    Becket Hui 2022/04
    
"""
import numpy as np
from miniball import get_bounding_ball
def calcCostFunc(ctrArr, volFctrArr):
    """
    calcCostFunc function computes the cost function of the given array of centroids,
    the cost function is based on the square distance and the volumetric factor
    input:
        ctrArr - array of centroids, structured as np[[x0, y0, z0],...]
        volFctrArr - array of volumetric factors as np[vf0, ...]
    return:
        cost - scalar cost
    """
    if ctrArr.shape[1] != 3:
        print('The input centroid array is not 3 dimensional.')
        return None
    if len(volFctrArr.shape) != 1:
        print('The input volume array is not 1 dimensional.')
        return None
    sqDist = calcSqDistFromIso(ctrArr)
    return np.sum(np.multiply(volFctrArr, sqDist))

def calcSqDistFromIso(ctrArr, isoCenter = None):
    """
    calcSqDistFromIso function computes the square distances between centroids and the isocenter
    input:
        ctrArr - array of centroids, structured as np[[x0, y0, z0],...]
        isoCenter - isoCenter, if no input, will calculate using calcIsoCenter
    return:
        sqDistArr - square distance for each centroids, structured as np[c0, ...]
    """
    if ctrArr.shape[1] != 3:
        print('The input array is not 3 dimensional.')
        return np.empty((0,3))
    if isoCenter is None:
        isoCenter = calcIsoCenter(ctrArr)
    distSqbyComp = np.square(ctrArr - isoCenter)
    return np.sum(distSqbyComp, axis=1)

def calcIsoCenter(ctrArr):
    """
    calcIsoCenter function computes isocenter from the array of given centroids using the minimum sphere method
    input:
        ctrArr - array of centroids, structured as np[[x0, y0, z0],...]
    return:
        isoArr - iso coordinate in the closest 0.5, structured as np[x, y, z]
    """
    if ctrArr.shape[1] != 3:
        print('The input array is not 3 dimensional.')
        return np.empty((0,3))
    cntr, r2 = get_bounding_ball(ctrArr)
    return 0.5*((2.0*cntr).round(decimals = 0))

# def calcIsoCenter(ctrArr):
#     """
#     calcIsoCenter function computes isocenter from the array of given centroids using the midpoint of the extreme method
#     input:
#         ctrArr - array of centroids, structured as np[[x0, y0, z0],...]
#     return:
#         isoArr - iso coordinate in the closest 0.5, structured as np[x, y, z]
#     """
#     if ctrArr.shape[1] != 3:
#         print('The input array is not 3 dimensional.')
#         return np.empty((0,3))
#     return 0.5*((ctrArr.min(axis = 0) + ctrArr.max(axis = 0)).round(decimals = 0))