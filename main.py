import numpy as np
from scipy.constants import *
import time
from qutip import *

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
N = 2                                                          #number of atoms

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

Dx = 0*ex                                           #Dipole moment vector
Dy = 0*ey
Dz = e*a0*ez                                        #Classical dipole moment in z-direction [C*m]
"""This might be the problem, as the transition dipole matrix element could be very different from the classical dipole moment."""
D = Qobj(Dz)                                        #Parallel polarization

#Vacuum decay rate for normalization of eigenenergies:
gamma_0 = (w0**3 * D.norm()**2) / (3 * pi * hbar * epsilon_0 * c**3)    #[Hz]
print(gamma_0)
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
    G = ((np.e**(k0*r_norm*complex(0,1)))/(4 * pi * k0**2 * r_norm**3)) * (
        (k0**2 * r_norm**2 + k0 * r_norm * complex(0, 1) - 1)*qeye(3) 
        + (-k0**2 * r_norm**2 - 3*k0*r_norm*complex(0, 1) + 3)*np.outer(r, r)/(r_norm**2)
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
    
    sigma_ge = Qobj([[0, 1], [0, 0]])                         #deexcitation
    sigma_eg = Qobj([[0, 0], [1, 0]])                         #excitation

    if i == 0:
        space = sigma_eg
    elif j == 0:
        space = sigma_ge
    else:
        space = qeye(2)
    for k in range(1, N):
        if k == i:
            space = tensor([space, sigma_eg])
        elif k == j:
            space = tensor([space, sigma_ge])
        else:
            space = tensor([space, qeye(2)])
        #print(k)

    return space
#Checking the values. They all seem reasonable, however.
#print(w0)          #384349305128205.1
#print(mu_0)        #1.2566370614359173e-06
#print(D.norm())    #1.602176634e-29
#print(gamma_0/w0)  #1.599244481542414e-10

def H_eff(N):
    H_eff = 0                                                       #Effective Hamiltonian of the system
    for i in range(N):                                              #Fill the effective Hamiltonian
        for j in range(N):
            if i == j:
                continue
            else:
                H_eff += (-mu_0 * w0**2) * D.trans() * Qobj(G[i,j]) * D * coherence_operators(i, j, N)
    return H_eff

def H(N, H_eff):
    sigma_ee = Qobj([[0, 0], [0, 1]])                             #excited state
    H = 0
    for i in range(N):

        if i == 0:
            space = sigma_ee
        else:
            space = qeye(2)
        for k in range(1, N):
            if k == i:
                space = tensor([space, sigma_ee])
            else:
                space = tensor([space, qeye(2)])

        H += hbar * w0 * space
        
    H += H_eff
    return H

"""
Section for producing plots
"""

Hamiltonian = H(N, H_eff(N))                #Hamiltionian as of eq. (5) in Asenjo-Garcia
energies = Hamiltonian.eigenenergies()
states = Hamiltonian.eigenstates()
#print(energies)

#Decay rates are the imaginary parts of the eigenenergies
decay_rates = - (2/hbar) * np.imag(energies) / gamma_0 #Normalized by vacuum decay rate
print(max(decay_rates))