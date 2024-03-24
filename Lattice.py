"""
Author: tzs820, s240670, Emil Henningsen

Compute position and interrelative displacement vectors of different lattice constructions
"""

from numpy import ndarray
from utils import *

class Lattice:

    #Declarations:

    positions = None
    displacements = None
    polarizations = None

    def __init__(self, pos : ndarray = None, disp : ndarray = None, pola : ndarray = None):
        """
        TODO: Description
        """
        self.positions = pos
        self.displacements = disp
        self.polarizations = pola

    """Construct different standard lattices:"""

    def linlat(self, N : int, d : float, direction : ndarray = None, polarizations : ndarray = None):
        """
        --- Parameters:

        N: integer, number of lattice sites
        d: float, distance between sites
        direction (optional): array (3) of floats, representing lattice vector.
            if none given, defaults to x direction.
        polarizations (optional): array (N x 3) of floats, representing each dipoles' polarization direction.
            if none given and internal container empty, defaults to z direction.

        --- Return:

        pos: array (N x 3) of floats, position vectors for each site
        rij: array (N x N x 3) of floats, interrelative displacement vectors
        """
        from numpy import zeros, linspace, array, full_like

        #if no direction is supplied, set standard direction to x:
        if direction == None:
            direction = array([1, 0, 0])
        
        #if no polarization is supplied and not already set, set standard to z.
        if polarizations == None and self.polarizations == None:
            polarizations = full_like(zeros(N), ez)    #orient all dipoles along z-direction.

        pos = zeros((N, 3))                                          #Array of position vectors of the lattice sites (atoms).
        ticks = linspace(-(N*d)/2, (N*d)/2, N)                       #displacement vectors in lattice
        for i in range(len(ticks)):
            pos[i, :] = direction * ticks[i]

        rij = zeros((N, N, 3))                                       #Array of interrelative displacement vectors between lattice sites.
        for i in range(N):                                              #Fill the array
            for j in range(i, N):
                rij[i, j] = pos[i] - pos[j]
                rij[j, i] = -rij[i, j]                               #this array is symmetric by construction.

        #Store in internal containers:
        self.positions = pos
        self.displacements = rij

    """Get methods:"""

    def getPositions(self) -> ndarray:
        return self.positions
    
    def getDisplacements(self) -> ndarray:
        return self.displacements
    
    def getPolarizations(self) -> ndarray:
        return self.polarizations
    
    """Set methods:"""
        
    def setPositions(self, pos : ndarray):
        self.positions = pos

    def setDisplacements(self, disp : ndarray):
        self.displacements = disp
    
    def setPolarizations(self, pola : ndarray):
        self.polarizations = pola