import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from qutip import *
from GreensTensor import *
import Hamiltonian
import Lattice
import Plots
from utils import *
from textwrap import wrap

#Program parameters:
N = 50                                                          #number of atoms


"""
Different inter-atomic distances in scalar (linear transverse) case.
"""

a = np.linspace(0.1, 2.0, 100)  #d/lambda0
d = 2*pi * a                    #The distance to feed G in units of 1/k0

lat = Lattice.Lattice()
scal = Hamiltonian.Hamiltonian()
pola = np.zeros((N,3))
for i in range(N):
    pola[i,:] = ex

p = Plots.Plots()
x = np.zeros(N)
for i in range(len(a)):
    x[:] = a[i]
    lat.linlat(N, d[i], ex, ex)
    displacements = lat.getDisplacements()
    G = fill_G(N, displacements)
    scal.block(N, G, pola)
    scal.eigenDecomposition()
    y = scal.getDecayRates()
    p.plot(x, y, 'b.')

plt.xlabel(r"$\frac{d}{\lambda_0}$", loc="right")
plt.ylabel(r"$\Gamma_\xi / \Gamma_0$", loc="top")
plt.ylim(0,3)
plt.title("\n".join(wrap(r"Varying $\frac{d}{\lambda_0} \in [0.1, 2.0]$ in linear chain of $N=$"+r"${}$".format(N)+" parallel polarized dipoles", 60)))
plt.savefig("figures/case_linear_parallel_var_distance_01_2.png", dpi=300)
p.show()