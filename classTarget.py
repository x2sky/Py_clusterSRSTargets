# -*- classTarget -*-#
"""
SYNOPSIS
    Class of a single target
    
DESCRIPTION
    Class that contains information of a single target
    
EXAMPLES
    target = SRSTarget(name, centroid coordinate, volume, volume factor)

VERSION 0.0
    first version
AUTHOR
    Becket Hui 2022/03
    
"""
import numpy as np
class SRSTarget:
    """
    SRS Target class that contains information of a SRS target
    """
    def __init__(self, name = '', cntr = np.zeros((3)), vol = 1.0, volFctr = 1.0):
        """
        SRSTarget class
        input (optional):
            name: target name
            cntr: centroid coordinate of the target in numpy array format
            vol: volume of target
            volFctr: volume factor
        """
        self.name = name
        self.centroid = cntr
        self.volume = vol
        self.volume_factor = volFctr