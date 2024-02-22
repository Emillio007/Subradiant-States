import numpy as np
from scipy.constants import *
import time
import matplotlib.pyplot as plt

def coherence_operators(i, j, N):
    
    if i == 0:
        space = np.array([[0, 0], [1, 0]])
    elif j == 0:
        space = np.array([[0, 1], [0, 0]])
    else:
        space = np.identity(2)
    for k in range(1, N):
        if k == i:
            space = np.kron(space, np.array([[0, 0], [1, 0]]))
        elif k == j:
            space = np.kron(space, np.array([[0, 1], [0, 0]]))
        else:
            space = np.kron(space, np.identity(2))
        #print(k)

    return space

N = 10

def doSum(N):
    sum = np.zeros((2**N, 2**N), dtype=float)
    for i in range(N):
        for j in range(N):
            sum += coherence_operators(i, j, N)
    return sum

def getNonZeroEntries(mat):
    return len(np.where(mat == 1)[0])

x = range(1, N+1)
y1 = [k**2 * k**2 for k in x]
y2 = [getNonZeroEntries(doSum(k)) for k in x]

plt.plot(x, y1, 'x', label="entries")
plt.plot(x, y2, 'x', label="non-zero entries")
plt.xscale('linear')
plt.yscale('linear')
plt.legend()
plt.show()

"""
Unfortunately, the scaling for the number of non-zero entries in the "sum of coherence operators" is also exponential.

In conclusion, I need to figure out how the "block diagonalization" works.
"""




















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









"""
Testing if scalar method and dyad tensor method gives equal results for linear lattice, dipoles polarized along z-direction.
"""
eigval, eigvec = np.linalg.eig(block)
eigval_scalar, eigvec_scalar = np.linalg.eig(G_scalar)
print(min(np.imag(eigval)))
print(min(np.imag(eigval_scalar)))
"""Scalar method and full dyad tensor method works equivalently! """