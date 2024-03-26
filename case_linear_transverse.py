import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from qutip import *
from GreensTensor import *
import Hamiltonian
import Lattice
import Plots
from utils import *

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
a = 100                                                         #d/lambda_0 = a
#d = 1/(2 * pi * a)                                              #dimensionless distance between the dipoles. The 2*pi might be wrong though????
d = 0.3                     #TODO: find ud af enheder for d og hvorfor størrelsen er, hvad den er. 
lattice = Lattice.Lattice()
lattice.linlat(N, d)        #initialize linear lattice
pos, rij, pola = lattice.getPositions(), lattice.getDisplacements(), lattice.getPolarizations()

p = Plots.Plots()
figDip = p.plotDipoles(lattice)

"""
Section for Hamiltonian:

All values are unitless, see notes from meeting 22/2. 
"""

#Scalar case:
scal = Hamiltonian.Hamiltonian()
scal.scalarham(N, rij)
scal.eigenDecomposition()

#Linear transverse case:
G = fill_G(N, rij)
block = Hamiltonian.Hamiltonian()
block.block(N, G, ez)       #initialize block hamiltonian with N dipoles and calculated G (vacuum) and ez pola direction.
block.eigenDecomposition()

#decay rates:
decay_rates = block.getDecayRates()
decay_rates_scalar = scal.getDecayRates()

"""Plotting decay rates of linear transverse and scalar case to see, if they are equal (which they should)"""
mytitle = r"$N = 50$ dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = 100$"
figDec = p.plotRatesLat(lattice, block, scalex="log", scaley = "log", title=None)
#figDecScal = p.plotRates(N, d, decay_rates_scalar, scaley="log", title="Scalar case")
#plt.savefig("figures/case_scalar_d_03.png", dpi=300)

p.show()

"""
Hermitiske del af Hamiltonian:
    (h + h^dagger)/2
"""

h = np.matrix(block.getHam())
herm = (h + h.getH())/2
herm_val, herm_vec = np.linalg.eigh(herm)
herm_diag = herm_vec.transpose() @ herm @ herm_vec
herm_trace = herm.trace()
herm_diag_trace = herm_diag.trace()

"""
Anti-hermitiske del af Hamiltonian:
    (h - h^dagger)/2
"""
anti = (h - h.getH())/2
anti_val, anti_vec = np.linalg.eig(anti)
anti_diag = np.matrix(anti_vec).getH() @ anti @ anti_vec
anti_trace = anti.trace()
anti_diag_trace = anti_diag.trace()

eigval_scalar, eigvec_scalar = block.getEigenDecomp()
print("Sum of eigenvalues: ", np.sum(eigval_scalar))
print("Trace of herm: ", herm_trace, " and trace of herm_diag: ", herm_diag_trace)
print("Trace of anti: ", anti_trace, " and trace of anti_diag: ", anti_diag_trace)
