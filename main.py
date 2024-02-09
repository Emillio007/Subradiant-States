import numpy as np
from scipy.constants import *
import time

#Declarations:
lambda_0 = 780e-9                                               #m (wavelength of the ??? transition)
w0 = 2*pi*c/lambda_0                                             #rad/s (angular frequency of the ??? transition)
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
        #rij[j, i] = -rij[i, j] no need to do this before G_0

"""
DIPOLE MOMENT
In this section, the dipole moment vector of chosen atom transition is calculated.
How does polarization direction of the dipole get into the calculation?
"""

#Dx = bra(psi0) @ q*ex @ ket(psi1)                                #Dipole moment vector of the atom transition
#Dy = bra(psi0) @ q*ey @ ket(psi1)
#Dz = bra(psi0) @ q*ez @ ket(psi1)
D = np.zeros(3)

"""
GREEN'S TENSOR (in free space)
In this section, the Green's tensor for the i'th dipole w.r.t. the j'th dipole is calculated.
The assumptions include:
 * Markovian regime: I.e. the bandwidth of atomic transition, w0, is small and the retardation of the field is negligible.
 * Free space: I.e. the Green's tensor is calculated in free space, where the dipoles are tightly trapped, i.e. d < lambda_0.
 * The stroung coupling regime is avoided: I.e. the electromagnetic field is not strongly coupled to the atomic transition (CQED).
 * The dipoles are point-like particles: I.e. r_i and r_j can be treated as fixed position vectors. 
"""

def G_0(r, w):
    """
    Green's tensor in free space
    """
    k0 = w/c
    r_norm = np.linalg.norm(r)
    G = ((np.e**(k0*r_norm*complex(0,1)))/(4*pi*epsilon_0*k0**2 * r_norm**3)) * (
        (k0**2 * r_norm**2 + k0*r_norm*complex(0, 1) - 1)*np.identity(3) 
        - (-k0**2 * r_norm**2 - 3*k0*r_norm*complex(0, 1) + 3)*np.outer(r, r)/(r_norm**2)
        )
    return G

G = np.zeros((N, N, 3, 3), dtype=complex)                        #Array of Green's tensors
for i in range(N):                                              #Fill the array of Green's tensors
    for j in range(i, N):
        if i == j:
            G[i, j] = np.zeros((3, 3), dtype=complex)
        else:
            G[i, j] = G_0(rij[i, j], w0)
            G[j, i] = G[i, j]

"""
EFFECTIVE HAMILTONIAN
In this section, the effective Hamiltonian of the system is calculated.
"""

def coherence_operators(i, j, N):
    """
    Coherence operators
    """
    space = np.identity(2)
    for k in range(N):
        if k == i:
            space = np.kron(space, np.array([[0, 0], [1, 0]]))
        elif k == j:
            space = np.kron(space, np.array([[0, 1], [0, 0]]))
        else:
            space = np.kron(space, np.identity(2))
        print(k)

time1 = time.time()
coherence_operators(8, 12, 20)
time2 = time.time()
print(time2-time1)

H_eff = np.zeros((N, N), dtype=complex)                          #Effective Hamiltonian of the system

for i in range(N):                                              #Fill the effective Hamiltonian
    for j in range(N):
        H_eff[i, j] += (-mu_0 * w0**2) * np.dot(D, np.dot(G[i, j], D))

print(H_eff.shape)