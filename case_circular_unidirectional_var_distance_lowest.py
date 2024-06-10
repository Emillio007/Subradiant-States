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
Different inter-atomic distances in circular unidir two-dim case.
"""
lowa, higha = 0.1, 0.4
a = np.linspace(lowa, higha, 100)  #d/lambda0
d = 2*pi * a                    #The distance to feed G in units of 1/k0

lowN, highN = 40, 60
N_var = range(lowN, highN+1)
a_static = 0.3
d_static = 2*pi * a_static   #d/lambda0 = 0.3

lat = Lattice.Lattice()
scal = Hamiltonian.Hamiltonian()
pola = np.zeros((N,3))
for i in range(N):
    pola[i,:] = ex

lowest_decay_rates = []
p = Plots.Plots()
def vary_d():
    for i in range(len(a)):
        lat.circlelat(N, d[i], distance_measure="inter", std_polarization="other", polarizations=pola)
        displacements = lat.getDisplacements()
        G = fill_G(N, displacements)
        scal.block(N, G, pola)
        scal.eigenDecomposition()
        y = scal.getDecayRates()
        lowest_decay_rates.append(y[0])

vary_d()
lowest_decay_rates_vary_d = np.array(lowest_decay_rates)

lowest_decay_rates_vary_N = []
def vary_N():
    for i in N_var:
        lat.circlelat(i, d_static, distance_measure="inter", std_polarization="other", polarizations=ex)
        displacements = lat.getDisplacements()
        pola = lat.getPolarizations()
        G = fill_G(i, displacements)
        scal.block(i, G, pola)
        scal.eigenDecomposition()
        y = scal.getDecayRates()
        lowest_decay_rates_vary_N.append(y[0])

vary_N()
lowest_decay_rates_vary_N = np.array(lowest_decay_rates_vary_N)

#third degree fits:
deg = 3
dega = 2
xs1 = np.linspace(lowa, higha, 1000)
xs2 = np.linspace(lowN, highN, 1000)
#fit1 = np.poly1d(np.polyfit(a, lowest_decay_rates_vary_d, dega))
interp1 = np.interp(xs1, a, lowest_decay_rates_vary_d)
#fit2 = scipy.optimize.curve_fit(lambda n: )
ys1 = interp1
#ys2 = fit2(xs2)

fig, ax1 = plt.subplots(layout="constrained")

#Plot first ax

ax2 = ax1.twiny()

ax1.plot(a, lowest_decay_rates_vary_d, 'x', label=r"var $\frac{d}{\lambda0}$")
#ax1.plot(xs1, ys1, 'r-', label="interpolation")
ax2.plot(N_var, lowest_decay_rates_vary_N, 'rx', label=r"var N")
#ax2.plot(xs2, ys2, '-', c="purple", label="3-deg fit")
ax2.grid(color="grey", linestyle="-.")
ax1.set_yscale("log")
ax1.set_xlabel(r"$\frac{d}{\lambda_0}$", loc="right")
ax2.set_xlabel("N", loc="right")
ax2.set_xticks([tick for tick in N_var if tick % 2 == 0])
ax1.set_ylabel(r"$\Gamma_{\xi=1} / \Gamma_0$", loc="top")
ax1.legend()
ax2.legend()

ax1.set_title("\n".join(wrap(r"Twin axes lowest decay rate vs. $\frac{d}{\lambda_0} \in [0.1, 0.4]$ for $N=%s$ and vs. $N \in [40, 60]$ for $\frac{d}{\lambda_0}=%s$ in circular chain of unidirectional polarized dipoles" % (N, a_static), 80)), loc="center")
plt.savefig("figures/case_circular_unidirectional_var_distance_01_04_var_N_40_60_lowest.png", dpi=300)
p.show()