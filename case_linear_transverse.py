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
a = 0.3                                                         #d/lambda_0 = a
#d = 1/(2 * pi * a)                                              #dimensionless distance between the dipoles. The 2*pi might be wrong though????
d = 100                     #TODO: find ud af enheder for d og hvorfor størrelsen er, hvad den er. 
lattice = Lattice()
lattice.linlat(N, d)        #initialize linear lattice
pos, rij, pola = lattice.getPositions(), lattice.getDisplacements(), lattice.getPolarizations()

p = Plots()
figDip, axDip = p.plotDipoles(lattice)

"""
Section for Hamiltonian:

All values are unitless, see notes from meeting 22/2. 
"""

#Scalar case:
G = fill_G(N, rij)
block = Hamiltonian()
block.block(N, G, ez)       #initialize block hamiltonian with N dipoles and calculated G (vacuum) and ez pola direction.

#decay rates:
decay_rates = block.getDecayRates()


#plt.savefig("figures/case_scalar.png", dpi=300)

plt.show()

"""
Hermitiske del af Hamiltonian:
    (h + h^dagger)/2
"""

h = np.matrix(block)
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

eigval_scalar = block.get
print("Sum of eigenvalues: ", np.sum(eigval_scalar))
print("Trace of herm: ", herm_trace, " and trace of herm_diag: ", herm_diag_trace)
print("Trace of anti: ", anti_trace, " and trace of anti_diag: ", anti_diag_trace)
