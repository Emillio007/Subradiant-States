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

p = Plots.Plots()

"""
Circular lattice, unidirectional polarization
"""

a = 0.1  #d/lambda0
d = 2*pi * a                    #The distance to feed G in units of 1/k0

pola = np.zeros((N,3))
for i in range(N):
    pola[i,:] = ex

lat = Lattice.Lattice()
lat.circlelat(N, d, distance_measure="inter", std_polarization="other", polarizations=pola)

mytitle = "\n".join(wrap(r"$N = %s$ dipoles in circular lattice with $\frac{d}{\lambda0}=%s$ polarized unidirectionally" % (N, a), 60))
p.plotDipolesPlane(lat, plane="xy", title=mytitle)

block = Hamiltonian.Hamiltonian()
pos, rij = lat.getPositions(), lat.getDisplacements()
G = fill_G(N, rij)
block.block(N, G, pola)

decay_rates = block.getDecayRates()
p.plotRatesLat(lat, block, scalex="log", scaley="log", title=mytitle)

p.show()
