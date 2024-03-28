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

angles = np.linspace(0, pi/2, 9)
#print(angles)
names = [r"0", r"$\pi/16$", r"$\pi/8$", r"$3\pi/16$", r"$\pi/4$", r"$5\pi/16$", r"$3\pi/8$", r"$7\pi/16$", r"$\pi/2$"]
fig, ax = plt.subplots(ncols=3, nrows=3)
fig.set_size_inches(20, 14)
indextable = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]
for angle, index, name in zip(angles, indextable, names):
    i, j = index
    dir = np.array([np.cos(angle), 0, np.sin(angle)])
    lat.linlat(N, d, ex, dir)
    displacements, pola = lat.getDisplacements(), lat.getPolarizations()
    G = fill_G(N, displacements)
    scal.block(N, G, pola)
    scal.eigenDecomposition()
    thistitle = r"Angle = %s" % name
    p.plotRatesLat(lat, scal, ax[i][j], scalex="log", scaley="log", title=thistitle)
    ax[i][j].set_xlim(19, 51)
    ax[i][j].set_ylim(0.1, 5)
    

fig.suptitle("\n".join(wrap(r"Linear chain of $N = %s$ dipoles with varying angle of polarization from x $ \in$"%N+" [%s, %s]" % (names[0], names[-1]), 120)))
plt.savefig("figures/case_linear_var_angle_0_pi05.png", dpi=300)
p.show()