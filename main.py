import numpy as np
from scipy.constants import *
import time

#Declarations:
lambda_0 = 780e-9                                               #m (wavelength of the ??? transition)
ex = np.array([1, 0, 0])                                         #x-axis unit vector
ey = np.array([0, 1, 0])                                         #y-axis unit vector
ez = np.array([0, 0, 1])                                         #z-axis unit vector

"""
Construct lattice
 (1) Finite linear chain
"""
d = 0.2 * lambda_0                                              #m (distance between the dipoles)
N = 20                                                          #number of atoms

pos = np.zeros((N, 3))                                          #Array of position vectors of the atoms
z_ticks = np.linspace(-(N*d)/2, (N*d)/2, N)                     #z-coordinates of the atoms
pos[:, 2] = z_ticks

rij = np.zeros((N, N, 3))                                       #Array of distance vectors between the atoms
for i in range(N):                                              #Fill the array of distance vectors
    for j in range(i, N):
        rij[i, j] = pos[i] - pos[j]
        rij[j, i] = -rij[i, j]

"""
DIPOLE MOMENT
In this section, the dipole moment vector of chosen atom transition is calculated.
How does polarization direction of the dipole get into the calculation?
"""

#Dx = bra(psi0) @ q*ex @ ket(psi1)                                #Dipole moment vector of the atom transition
#Dy = bra(psi0) @ q*ey @ ket(psi1)
#Dz = bra(psi0) @ q*ez @ ket(psi1)

"""
GREEN'S TENSOR (in free space)
In this section, the Green's tensor for the i'th dipole w.r.t. the j'th dipole is calculated.
The assumptions include:
 * Markovian regime: I.e. the bandwidth of atomic transition, w0, is small and the retardation of the field is negligible.
 * Free space: I.e. the Green's tensor is calculated in free space, where the dipoles are tightly trapped, i.e. d < lambda_0.
 * The stroung coupling regime is avoided: I.e. the electromagnetic field is not strongly coupled to the atomic transition (CQED).
 * The dipoles are point-like particles: I.e. r_i and r_j can be treated as fixed position vectors. 
"""

