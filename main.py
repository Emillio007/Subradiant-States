import numpy as np
from scipy.constants import *
import time
from qutip import *
from GreensTensor import *
from Hamiltonian import *

#Declarations:
lambda_0 = 780e-9                                               #[m] (wavelength of the ??? transition)
w0 = c/lambda_0                                             #[1/s] (frequency of the ??? transition)
ex = np.array([1, 0, 0])                                         #x-axis unit vector
ey = np.array([0, 1, 0])                                         #y-axis unit vector
ez = np.array([0, 0, 1])                                         #z-axis unit vector
a0 = 1e-10                                                     #[m] (atomic unit of length)


#N.B.:
#hbar = 1


"""
Construct lattice
 (1) Finite linear chain along z-axis
"""
d = 0.2 * lambda_0                                              #m (distance between the dipoles)
N = 10                                                          #number of atoms

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

Dx = 0*ex                                           #Dipole moment vector
Dy = 0*ey
Dz = e*a0*ez                                        #Classical dipole moment in z-direction [C*m]
"""This might be the problem, as the transition dipole matrix element could be very different from the classical dipole moment."""
D = Qobj(Dx+Dy+Dz)                                        #Parallel polarization
Dnorm = D.norm()
print(Dnorm)
#Vacuum decay rate for normalization of eigenenergies:
gamma_0 = (w0**3 * Dnorm**2) / (3 * pi * hbar * epsilon_0 * c**3)    #[Hz]

G = fill_G(N, rij, w0, type="free_space")

"""
Section for Hamiltonian:
"""

#We want Hamiltonian in basis of {|e_j>}, meaning one excitation on j'th subspace.
#Do this by computing every matrix element of NxN matrix. It is the block of one excitation in full Hamiltonian
block = np.zeros((N,N))
for i in range(N):
    for j in range(N):
        if i == j:
            H_eff = 1
        else:
            H_eff = -mu_0 * w0**2 * ez.transpose() @ G[i,j] @ ez
        block.data[i,j] += H_eff
    block.data[i,i] += hbar * w0

print(block)
"""
Section for producing plots:
"""

#Hamiltonian = H(N, H_eff(N))                #Hamiltionian as of eq. (5) in Asenjo-Garcia
#energies = Hamiltonian.eigenenergies()
#states = Hamiltonian.eigenstates()
#print(energies)

#Decay rates are the imaginary parts of the eigenenergies
#decay_rates = - (2/hbar) * np.imag(energies) / gamma_0 #Normalized by vacuum decay rate
#print(Hamiltonian)


"""
Testing out equation 4.b'
"""

#normalized decay rates: gamma_ij / gamma_0
gamma = np.zeros((N, N))
fac = (6 * np.pi * epsilon_0 * mu_0 * c**3) / (w0)
for i in range(N):
    for j in range(N):
        gamma[i,j] += fac * ez.transpose() @ np.imag(G[i,j]) @ ez   #Using z direction vector for simplicity.

print(np.min(gamma))