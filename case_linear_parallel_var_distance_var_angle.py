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
Different inter-atomic distances in linear var angle case.
"""
ncases = 10
a = np.linspace(0.2, 0.4, ncases)  #d/lambda0
d = 2*pi * a                    #The distance to feed G in units of 1/k0

#Unit vector with angle for polarization (respective to linear chain axis, ex)
angle = 5*pi/16
pola_vec = np.array([np.cos(angle), 0, np.sin(angle)])

lat = Lattice.Lattice()
scal = Hamiltonian.Hamiltonian()
pola = np.zeros((N,3))
for i in range(N):
    pola[i,:] = pola_vec

p = Plots.Plots()
x = np.zeros(N)

#x axis for the last 30 decay rate points
xax = range(N)[-30:-1]

#Set color cycle to gradient:
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.gray(np.linspace(0, 1, ncases)))

for i in range(len(a)):
    x[:] = a[i]
    lat.linlat(N, d[i], ex, pola)
    displacements = lat.getDisplacements()
    G = fill_G(N, displacements)
    scal.block(N, G, pola)
    scal.eigenDecomposition()
    y = scal.getDecayRates()
    yax = y[-30:-1]
    p.plot(xax, yax, '.')
    

plt.xlabel(r"$\xi$", loc="right")
plt.ylabel(r"$\Gamma_\xi / \Gamma_0$", loc="top")
plt.title("\n".join(wrap(r"Varying $\frac{d}{\lambda_0} \in [%s, %s]$ in linear chain of $N=%s$ dipoles polarized w $\theta=%s$" % (a[0], a[-1], N, angle), 80)), loc="center")
#plt.savefig("figures/case_linear_parallel_var_distance_01_2.png", dpi=300)
p.show()