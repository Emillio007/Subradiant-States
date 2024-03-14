import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
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

x = pos[:,0]
z = np.full_like(x, 0)
polax = np.full_like(x, 0)  #polarization direction on plot, x
polaz = np.full_like(x, 1)  #-||-, z
plt.figure()
plt.plot(x, z, 'o', color="black", label="sites")
plt.quiver(x, z, polax, polaz, scale=15, width=0.005, color="red", label=r"$\hat{d}$", pivot="mid")
plt.ylim(-1, 1)
plt.xlabel(r"$\mathbf{\hat{x}}$", loc="right")
plt.ylabel(r"$\mathbf{\hat{z}}$", loc="top")
plt.title(r"Linear lattice of $N=50$ dipoles, $\frac{d}{\lambda_0}=0.3$")
plt.legend()

G = fill_G(N, rij)              

"""
Section for Hamiltonian:

All values are unitless, see notes from meeting 22/2. 
"""

#Scalar case:
block = block(N, G, ez)

#Eigenvalues:
eigval_scalar, eigvec_scalar = np.linalg.eig(block)
decay_rates = 2 * np.imag(eigval_scalar)                       #Factor of 2, see Asenjo-Garcia et al.
decay_rates.sort()

plt.figure()
plt.plot(range(1, N+1), decay_rates, 'o')
plt.yscale("linear")
plt.xscale("linear")
plt.xlabel(r"$\mathbf{\xi \in [1,N]}$", loc="right")
plt.ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
plt.title(r"N=50 dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = 0.3$")
#plt.savefig("figures/case_scalar.png", dpi=300)


plt.show()

"""
Hermitiske del af Hamiltonian:
    (h + h^dagger)/2
"""

h = np.matrix(block)
herm = (h + h.getH())/2
herm_val, herm_vec = np.linalg.eigh(herm)
herm_diag = herm_vec.transpose() @ herm @ herm_vec
herm_trace = herm.trace()
herm_diag_trace = herm_diag.trace()

"""
Anti-hermitiske del af Hamiltonian:
    (h - h^dagger)/2
"""
anti = (h - h.getH())/2
anti_val, anti_vec = np.linalg.eig(anti)
anti_diag = np.matrix(anti_vec).getH() @ anti @ anti_vec
anti_trace = anti.trace()
anti_diag_trace = anti_diag.trace()

print("Sum of eigenvalues: ", np.sum(eigval_scalar))
print("Trace of herm: ", herm_trace, " and trace of herm_diag: ", herm_diag_trace)
print("Trace of anti: ", anti_trace, " and trace of anti_diag: ", anti_diag_trace)
