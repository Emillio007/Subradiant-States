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
N = 50                                                          #number of atoms

#Declarations:


"""
Construct lattice
 (1) Finite linear chain along x-axis
"""
a = 1/100000                          #d/lambda_0 = a
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

#decay rates:
decay_rates = block.getDecayRates()

mytitle = "\n".join(wrap(r"$N = $" + "{}".format(N) + r" dipoles in linear lattice, polarized in x-direction, $\frac{d}{\lambda_0} = $" + f"{a}", 60))

figDip, axDip = p.plotDipolesPlane(lattice, plane = "xz", title=mytitle, xlim=None, ylim=(-1,1), ham=block, index=block.getSortedIndex()[0])
figDip.set_size_inches(12, 5)
#plt.savefig("figures/dipoles_case_linear_parallel_d_03.png", dpi=300)

eigval, eigvec = block.getEigenDecomp()

def check_expectation_value(index : int) -> None:
    ind = block.getSortedIndex()[index]
    vec = eigvec[:,ind]
    #print(vec)
    expec_val = np.transpose(np.conjugate(vec)) @ block.getHam() @ vec
    print("Actual eigenvalue: ", eigval[ind], "\n Expectation value: ", expec_val)

#check_expectation_value(0)
#Checking the expectation value against the eigenvalue returns the same. So still positive imaginary component for sufficient conditions.

"""Plotting decay rates of linear parallel """

figDec = p.plotRatesLat(lattice, block, scalex="log", scaley = "log", title=mytitle)
#plt.savefig("figures/case_linear_parallel_d_03.png", dpi=300)
#print(decay_rates)
#p.show()
