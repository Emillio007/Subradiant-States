"""
Author: tzs820, s240670, Emil Henningsen

Compute position and interrelative displacement vectors of different lattice constructions
"""

from numpy import ndarray
from typing import Literal
from warnings import warn

class Lattice:

    #Declarations:
    
    N = None
    d = None

    positions = None
    displacements = None
    polarizations = None

    #Warnings and exceptions:
    warning_pola_none = Warning("No polarization is provided, proceeds with ez for all dipoles. ")
    exception_not_yet_implemented = Exception("The requested feature has not been implemented yet. ")

    def __init__(self, N : int = None, d : float = None, pos : ndarray = None, disp : ndarray = None, pola : ndarray = None):
        """
        TODO: Description
        """
        self.positions = pos
        self.displacements = disp
        self.polarizations = pola
        self.N = N
        self.d = d

    def clearLat(self) -> None:
        self.positions = None
        self.displacements = None
        self.polarizations = None
        self.N = None
        self.d = None

    """Get methods:"""

    def getN(self) -> int:
        return self.N
    
    def getd(self) -> float:
        return self.d

    def getPositions(self) -> ndarray:
        return self.positions
    
    def getDisplacements(self) -> ndarray:
        return self.displacements
    
    def getPolarizations(self) -> ndarray:
        return self.polarizations
    
    """Set methods:"""

    def setN(self, N : int) -> None:
        self.N = N
    
    def setd(self, d : float)-> None:
        self.d = d
        
    def setPositions(self, pos : ndarray)-> None:
        self.positions = pos

    def setDisplacements(self, disp : ndarray)-> None:
        self.displacements = disp
    
    def setPolarizations(self, pola : ndarray)-> None:
        self.polarizations = pola

    """Helper functions:"""

    def fillPola(self, N, dir : ndarray) -> ndarray:
        """
        TODO: Description
        """
        from numpy import zeros
        polarizations = zeros((N,3))
        for i in range(N):
            polarizations[i,:] = dir
        return polarizations
    
    def fillDisplacements(self, N : int, pos : ndarray) -> ndarray:
        """
        TODO: Description
        """
        from numpy import zeros

        rij = zeros((N, N, 3))                                       #Array of interrelative displacement vectors between lattice sites.
        for i in range(N):                                              #Fill the array
            for j in range(i, N):
                rij[i, j] = pos[i] - pos[j]
                rij[j, i] = -rij[i, j]                               #this array is symmetric by construction.

        return rij

    """Construct different standard lattices:"""

    def linlat(self, N : int, d : float, direction : ndarray = None, polarizations : ndarray = None)-> None:
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
        from utils import ex, ez

        #If already init, re-init.
        self.clearLat()

        #if no direction is supplied, set standard direction to x:
        if direction is None:
            direction = ex
        
        #if no polarization is supplied and not already set, set standard to z.
        if polarizations is None and self.polarizations is None:
            polarizations = self.fillPola(N, ez)                #Orient all dipoles along z
        elif polarizations.shape == (3,):
            dir = polarizations
            polarizations = self.fillPola(N, dir)

        pos = zeros((N, 3))                                          #Array of position vectors of the lattice sites (atoms).
        ticks = linspace(-(N*d)/2, (N*d)/2, N)                       #displacement vectors in lattice
        for i in range(len(ticks)):
            pos[i, :] = direction * ticks[i]

        rij = self.fillDisplacements(N, pos)

        #Store in internal containers:
        self.setN(N)
        self.setd(d)
        self.setPositions(pos)
        self.setDisplacements(rij)
        self.setPolarizations(polarizations)

    def circlelat(self, N : int, d : float, distance_measure : Literal["inter", "radius"] = "inter", 
                  std_polarization : Literal["inwards", "inwards alternating", "outwards", "outwards alternating", "azimuthal",
                                             "azimuthal alternating", "other"] = "other", polarizations : ndarray = None) -> None:
        """
        TODO: Description

        if std_polarization == "other", polarizations has to be supplied (either a single vec (3,) or full array (N,3)).
        """
        from numpy import zeros, pi, cos, sin, sqrt
        from utils import ez

        #If already init, clear lattice:
        self.clearLat()

        #do cirlce:
        def fillCircle(N : int, r : float, angle : float) -> ndarray:
            """
            TODO: Description
            """
            pos = zeros((N,3))
            for i in range(N):
                #fill positions for dipoles in xy-plane
                pos[i,0] = r * cos(i * angle)
                pos[i,1] = r * sin(i * angle)
            return pos
        b = 0                                           #inter-dipole distance
        r = 0                                           #radius
        angle = 2*pi / N
        match(distance_measure):
            case "inter":
                b = d
                r = b / sqrt(2 * (1 - cos(angle)))      #Cosine relation
            case "radius":
                r = d
                b = r * sqrt(2 * (1 - cos(angle)))      #Cosine relation
        pos = fillCircle(N, r, angle)

        #Displacements:
        rij = self.fillDisplacements(N, pos)

        #Do polarization:
        pola = None
        match(std_polarization):
            case "other":
                if polarizations is None and self.polarizations is None:
                    #Raise warning and choose standard direction
                    warn(self.warning_pola_none)
                    pola = self.fillPola(N, ez)
                elif polarizations.shape == (3,):
                    #provided polarization is a single vector, fill entire array.
                    pola = self.fillPola(N, polarizations)
                else:
                    #Last case, the entire array is supplied:
                    pola = polarizations
            case "inwards":
                #some vector math..
                pass
            case "inwards alternating":
                raise self.exception_not_yet_implemented
            case "outwards":
                raise self.exception_not_yet_implemented
            case "outwards alternating":
                raise self.exception_not_yet_implemented
            case "azimuthal":
                raise self.exception_not_yet_implemented
            case "azimuthal alternating":
                raise self.exception_not_yet_implemented

        #Store in internal containers:
        self.setN(N)
        self.setd(d)
        self.setPositions(pos)
        self.setDisplacements(rij)
        self.setPolarizations(pola)