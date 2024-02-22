"""
Author: tzs820, s240670, Emil Henningsen

Compute position and interrelative displacement vectors of different lattice constructions
"""

def linlat(N, d, direction):
    """
    
    """
    from numpy import zeros, linspace

    pos = zeros((N, 3))                                          #Array of position vectors of the lattice sites (atoms).
    ticks = linspace(-(N*d)/2, (N*d)/2, N)                       #displacement vectors in lattice
    for i in range(len(ticks)):
        pos[i, :] = direction * ticks[i]

    rij = zeros((N, N, 3))                                       #Array of interrelative displacement vectors between lattice sites.
    for i in range(N):                                              #Fill the array
        for j in range(i, N):
            rij[i, j] = pos[i] - pos[j]
            rij[j, i] = -rij[i, j]                               #this array is symmetric by construction.

    return pos, rij