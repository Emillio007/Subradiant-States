import numpy as np
from numpy import min
import matplotlib.pyplot as plt
from scipy.constants import *
from qutip import *
from GreensTensor import *
import Hamiltonian
import Lattice
import Plots
from utils import *
from textwrap import wrap

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

#Program parameters:
N = 50                                                          #number of atoms
npoints = 20                                                    #number of data points

"""
Construct lattice
 (1) Finite linear chain along x-axis
"""
a = 0.3                          #d/lambda_0 = a
d = 2*pi * a                     #Faktor 2pi fordi r i enheder af 1/k0
theta = np.pi/2
lattice = Lattice.Lattice()
block = Hamiltonian.Hamiltonian()

def func(angle):
    "Return minimum decay rate by angle"

    lattice.twopiece(N, d, angle, ex, ex)        #initialize linear lattice
    pos, rij, pola = lattice.getPositions(), lattice.getDisplacements(), lattice.getPolarizations()
        
    #Linear parallel case:
    G = fill_G(N, rij)
    block.block(N, G, pola)       #initialize block hamiltonian with N dipoles and calculated G (vacuum) and ez pola direction.
    block.eigenDecomposition()

    #decay rates:
    decay_rates = block.getDecayRates()

    return min(decay_rates)

angles = np.linspace(0, theta, npoints)
rates = []
for angle in angles:
    rates.append(func(angle))
rates = np.array(rates)

p = Plots.Plots()

mytitle = "\n".join(wrap(r"$N = $" + "{}".format(N) + r" dipoles in broken linear lattice, polarized in x-direction, $\frac{d}{\lambda_0} = $" + f"{a}", 60))

#plot:
plt.plot(angles, rates, 'o')
plt.yscale("log")
plt.show()

#figDip, axDip = p.plotDipolesPlane(lattice, plane = "xz", title=mytitle, xlim=None, ylim=None, ham=block, index=block.getSortedIndex()[0])
#figDip.set_size_inches(12, 5)
#plt.savefig("figures/dipoles_case_linear_parallel_d_03.png", dpi=300)

"""Plotting decay rates of linear parallel """
#print(decay_rates)
#figDec = p.plotRatesLat(lattice, block, scalex="log", scaley = "log", title=mytitle)
#plt.savefig("figures/case_linear_parallel_d_03.png", dpi=300)

#p.show()
