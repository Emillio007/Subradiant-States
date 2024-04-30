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

a = 0.3  #d/lambda0 .... Bug: for a = 0.01 the lowest decay rate is negative.
d = 2*pi * a                    #The distance to feed G in units of 1/k0

pola = np.zeros((N,3))
for i in range(N):
    pola[i,:] = ex

lat = Lattice.Lattice()
lat.circlelat(N, d, distance_measure="inter", std_polarization="other", polarizations=pola)

block = Hamiltonian.Hamiltonian()
pos, rij = lat.getPositions(), lat.getDisplacements()
G = fill_G(N, rij)
block.block(N, G, pola)

decay_rates = block.getDecayRates()
print(decay_rates)

mytitle = "\n".join(wrap(r"$N = %s$ dipoles in circular lattice with $\frac{d}{\lambda0}=%s$ polarized unidirectionally" % (N, a), 60))
figDip, axDip = p.plotDipolesPlane(lat, plane="xy", title=mytitle, ham=block, index=block.getSortedIndex()[0])

p.plotRatesLat(lat, block, scalex="log", scaley="log", title=mytitle)

p.show()
