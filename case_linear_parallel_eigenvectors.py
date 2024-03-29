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

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

#Program parameters:
N = 20                                                          #number of atoms

#Declarations:


"""
Construct lattice
 (1) Finite linear chain along x-axis
"""
a = 0.3                          #d/lambda_0 = a
d = 2*pi * a                     #Faktor 2pi fordi r i enheder af 1/k0
lattice = Lattice.Lattice()
lattice.linlat(N, d, ex, ex)        #initialize linear lattice
pos, rij, pola = lattice.getPositions(), lattice.getDisplacements(), lattice.getPolarizations()
p = Plots.Plots()

"""
Section for Hamiltonian:

All values are unitless, see notes from meeting 22/2. 
"""

#Linear parallel case:
G = fill_G(N, rij)
block = Hamiltonian.Hamiltonian()
block.block(N, G, pola)       #initialize block hamiltonian with N dipoles and calculated G (vacuum) and ez pola direction.
block.eigenDecomposition()
val, vec = block.getEigenDecomp()

#decay rates:
decay_rates = block.getDecayRates()

figDip, axDip = None, None
amount = 10         #Number of lowest eigenstates to plot
mytitle = "\n".join(wrap(r"The %s highest first-excitation eigenvectors of $N = %s$ dipoles in linear lattice, polarized in x-direction, $\frac{d}{\lambda_0} = %s$" % (amount, N, a), 120))
positions = pos #Original positions
startind = N - amount
for i in range(amount):
    positions[:,2] = startind + i                 #Add i to all z-coordinates to shift the i'th dipole chain upwards.
    lattice.setPositions(positions)    
    figDip, axDip = p.plotDipolesPlane(lattice, ax=axDip, plane = "xz", title=mytitle, xlim=None, ylim=None, ham=block, index=block.getSortedIndex()[startind+i], legend=False)
    print("Decay rate of ", startind+i, "'th state: ", 2*np.imag(val[startind+i]))
    """
    TODO: Find ud af, hvorfor argsort ikke sorterer rigtigt. Altså, hvorfor er de sidste fire tilstande sorteret på hovedet? 
    """
print(decay_rates)
axDip.set_ylabel(r"$|\psi_i> \in [%s, %s)$"%(startind, N), loc="top")
figDip.set_size_inches(14, 5)
#plt.savefig("figures/dipoles_case_linear_parallel_10higheststates.png", dpi=300)

"""Plotting decay rates of linear parallel """

#figDec = p.plotRatesLat(lattice, block, scalex="log", scaley = "log", title=mytitle)
#plt.savefig("figures/case_linear_parallel_d_03.png", dpi=300)

p.show()
