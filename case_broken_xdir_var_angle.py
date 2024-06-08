import numpy as np
from numpy import min
import matplotlib.pyplot as plt
import matplotlib
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
npoints = 50                                                    #number of data points

"""
Construct lattice
 (1) Finite linear chain along x-axis
"""
a = 0.3                          #d/lambda_0 = a
d = 2*pi * a                     #Faktor 2pi fordi r i enheder af 1/k0
theta = np.pi/2
lattice = Lattice.Lattice()
block = Hamiltonian.Hamiltonian()
p = Plots.Plots()

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

mytitle = "\n".join(wrap(r"$N = %s$ dipoles in broken linear lattice, polarized in x-direction, $\frac{d}{\lambda_0} = %s$" % (N, a), 80))

fig, ax = plt.subplots(nrows = 3)
angles = np.linspace(0, theta, npoints)
rates = []
index = 0
pi = 0
amplRange = np.array([])
for angle in angles:
    rates.append(func(angle))
    if index == 0 or index == 24 or index == 49:
        f, a, ampl = p.plotDipolesPlane(lattice, ax=ax[pi], plane="xz", title=mytitle, xlim=None, ylim=None, ham=block, index=block.getSortedIndex()[0], drawColorbar = False, returnAmpls=True)
        pi += 1
        amplRange = np.append(amplRange, ampl)
    index += 1
rates = np.array(rates)
print(rates)

fig.subplots_adjust(right=0.8)
cmap = matplotlib.cm.coolwarm
norm = matplotlib.colors.Normalize(vmin=min(amplRange), vmax=max(amplRange))
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
fig.colorbar(sm, cax=fig.add_axes([0.85, 0.15, 0.05, 0.7]), label=r"$|c_j|$")
#fig.set_layout_engine("constrained")
#plt.savefig("figures/dipoles_case_linear_broken.png", dpi=300)

#plot:
fig = plt.figure()
fig.set_size_inches(10, 6)
mytitle = "\n".join(wrap(r"Decay rates of most subradiant modes of $N = %s$ dipoles in broken (with angle $\theta$) linear lattice" % (N), 60))
plt.plot(angles, rates, 'o', label="Decay rates")
plt.yscale("log")
plt.title(mytitle)
plt.xlabel(r"$\theta \in [0, \frac{\pi}{2}]$", loc="right")
plt.ylabel(r"$\Gamma_{\xi=1} / \Gamma_0$", loc="top")
plt.legend()
#plt.savefig("figures/case_linear_broken_rates.png", dpi=300)
plt.show()

#figDip, axDip = p.plotDipolesPlane(lattice, plane = "xz", title=mytitle, xlim=None, ylim=None, ham=block, index=block.getSortedIndex()[0])
#figDip.set_size_inches(12, 5)
#plt.savefig("figures/dipoles_case_linear_parallel_d_03.png", dpi=300)

"""Plotting decay rates of linear parallel """
#print(decay_rates)
#figDec = p.plotRatesLat(lattice, block, scalex="log", scaley = "log", title=mytitle)
#plt.savefig("figures/case_linear_parallel_d_03.png", dpi=300)

#p.show()
