import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
import time
from qutip import *
from GreensTensor import *
from Hamiltonian import *
from Lattice import *

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

#Program parameters:
N = 50                                                          #number of atoms

#Declarations:
ex = np.array([1, 0, 0])                                         #x-axis unit vector
ey = np.array([0, 1, 0])                                         #y-axis unit vector
ez = np.array([0, 0, 1])                                         #z-axis unit vector

"""
Construct lattice
 (1) Finite linear chain along x-axis
"""
a = 0.3                                                         #d/lambda_0 = a
d = 1/(2 * pi * a)                                              #dimensionless distance between the dipoles. The 2*pi might be wrong though????

pos, rij = linlat(N, d, ex)

G = fill_G(N, rij)

"""
Section for Hamiltonian:

All values are unitless, see notes from meeting 22/2. 
"""

#We want Hamiltonian in basis of {|e_j>}, meaning one excitation on j'th subspace.
#Do this by computing every matrix element of NxN matrix. It is the block of one excitation in full Hamiltonian
block = np.zeros((N,N), dtype=complex)

for i in range(N):
    for j in range(N):
        if i == j:
            block[i, j] += 1
        else:
            block[i, j] += ez.transpose() @ G[i,j] @ ez        #Dipoles polarized along z-direction.  

#Scalar case:
G_scalar = scalar(N, rij)

#Eigenvalues:
eigval_scalar, eigvec_scalar = np.linalg.eig(G_scalar)
decay_rates = 2 * np.imag(eigval_scalar)                       #Factor of 2, see Asenjo-Garcia et al.
decay_rates.sort()

plt.plot(range(1, N+1), decay_rates, 'o')
plt.yscale("log")
plt.xlabel(r"$\mathbf{\xi \in [1,N]}$", loc="right")
plt.ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
plt.show()