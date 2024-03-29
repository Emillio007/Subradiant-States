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

a = 0.3  #d/lambda0
d = 2*pi * a                    #The distance to feed G in units of 1/k0

lat = Lattice.Lattice()
lat.circlelat(N, d, distance_measure="inter", std_polarization="inwards", polarizations=None)

block = Hamiltonian.Hamiltonian()
pos, rij, pola = lat.getPositions(), lat.getDisplacements(), lat.getPolarizations()
G = fill_G(N, rij)
block.block(N, G, pola)

decay_rates = block.getDecayRates()
print(decay_rates)

#Plotting dipoles in plane with colormap of probability amplitude norms of most subradiant eigenstate
mytitle = "\n".join(wrap(r"$N = %s$ dipoles in circular lattice with $\frac{d}{\lambda0}=%s$ polarized inwards" % (N, a), 60))
p.plotDipolesPlane(lat, plane="xy", title=mytitle, ham=block, index=block.getSortedIndex()[0])  #Indexes have been sorted, bc getDecayRates have been called above.

p.plotRatesLat(lat, block, scalex="log", scaley="log", title=mytitle)

p.show()
