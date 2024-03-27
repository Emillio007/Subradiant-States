import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from qutip import *
from GreensTensor import *
import Hamiltonian
import Lattice
import Plots
from utils import *

#Program parameters:
N = 50                                                          #number of atoms


"""
Different inter-atomic distances in scalar (linear transverse) case.
"""

d = np.linspace(0.1, 2.0, 100)

lat = Lattice.Lattice()
scal = Hamiltonian.Hamiltonian()

p = Plots.Plots()
x = np.zeros(N)
for distance in d:
    x[:] = distance
    lat.linlat(N, distance)
    displacements = lat.getDisplacements()
    G = fill_G(N, displacements)
    scal.block(N, G, ez)
    scal.eigenDecomposition()
    y = scal.getDecayRates()
    p.plot(x, y, 'o', color="blue")

plt.xlabel("d")
plt.ylabel(r"$\Gamma_\xi / \Gamma_0$")
plt.title(r"Varying $d \in [0.1, 2.0]$")
#plt.savefig("figures/case_scalar_var_distance_10_200.png", dpi=300)
p.show()