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
Different polarization angles distances in linear transverse (scalar) case.
"""

a = 0.3
d = 2*pi * a

lat = Lattice.Lattice()
scal = Hamiltonian.Hamiltonian()

p = Plots.Plots()

angles = np.linspace(0, pi/2, 6)
fig, ax = plt.subplots(ncols=3, nrows=2)
fig.set_size_inches(20, 10)
indextable = [(0,0), (1,0), (2,0), (0,1), (1,1), (2,1)]
for angle, index in zip(angles, indextable):
    j, i = index
    dir = np.array([np.cos(angle), 0, np.sin(angle)])
    lat.linlat(N, d, ex, dir)
    displacements, pola = lat.getDisplacements(), lat.getPolarizations()
    G = fill_G(N, displacements)
    scal.block(N, G, pola)
    scal.eigenDecomposition()
    thistitle = r"Angle = %s" % angle
    p.plotRatesLat(lat, scal, ax[i][j], scalex="log", scaley="log", title=thistitle)
    

fig.suptitle("\n".join(wrap(r"Linear chain of $N = %s$ dipoles with varying angle of polarization from x $ \in [%s, %s]$" % (N, angles[0], angles[-1]), 120)))
#plt.savefig("figures/case_scalar_var_distance_10_200.png", dpi=300)
p.show()