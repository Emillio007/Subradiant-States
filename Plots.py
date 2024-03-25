"""
Author: tzs820, s240670, Emil Henningsen

Plotting module
"""

import matplotlib.pyplot as plt
import Lattice
import Hamiltonian
from numpy import ndarray
from cycler import cycler

class Plots:

    figures = {}

    def __init__(self, figures : dict = None):
        """
        Initialize plot manager. 
        """

        if not figures == None:
            self.figures = figures

        #Upon init, set matplotlib standards. 
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica",
            "axes.labelsize": 16,
            "axes.labelweight": "bold",
            "axes.titlesize": 16,
            "axes.titlelocation": "left",
            "axes.prop_cycle": cycler(color=["black", "red", "blue", "purple"]),
            "axes.grid": True
        })

    """modules:"""

    def add(self, fig, name : str = None):
        """
        Add a figure to the instantiated manager container with given name (dict key). 
        """
        if name == None:
            number = len(self.figures)
            name = number
        self.figures[name] = fig

    def remove(self, name : str):
        """
        Remove specified figure from container (name = dict key)
        """
        self.figures.pop(name)

    def show(self):
        plt.show()

    """SET modules:"""

    def setInteractive(self, interactive : bool = False):
        """
        Set matplotlib mode to interactive. If false (standard), figures are not shown until show() is called.
        If true, figures are shown immediately after creation.
        """
        
        if interactive == False:
            plt.ioff()
        else:
            plt.ion()

    """GET modules:"""

    def isInteractive() -> bool:
        return plt.isinteractive()
            
    """Different standard type plots: """
    def plotDipoles(self, lat : Lattice) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description
        """
        
        pos, dir, pola = lat.getPositions(), lat.getDisplacements(), lat.getPolarizations()

        x = pos[:,0]    #All x coordinates
        y = pos[:,1]    #All y coordinates
        z = pos[:,2]    #All z coordinates

        polax = pola[:,0]
        polay = pola[:,1]
        polaz = pola[:,2]

        fig, ax = plt.figure()
        ax.plot(x, z, 'o', color="black", label="sites")
        ax.quiver(x, z, polax, polaz, scale=15, width=0.005, color="red", label=r"$\hat{d}$", pivot="mid")
        ax.ylim(-1, 1)
        ax.xlabel(r"$\mathbf{\hat{x}}$", loc="right")
        ax.ylabel(r"$\mathbf{\hat{z}}$", loc="top")
        ax.title(r"Linear lattice of $N=50$ dipoles, $\frac{d}{\lambda_0}=0.3$")
        ax.legend()
        
        return fig, ax
    
    #Plot rates manually
    def plotRates(self, N : int, rates : ndarray) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description

        Find smart way to transfor information from hamiltonian to plot. Maybe actually just feed hamiltonian directly?
        """
        fig, ax = plt.figure()
        ax.plot(range(1, N+1), rates, 'o')
        ax.yscale("linear")
        ax.xscale("linear")
        ax.xlabel(r"$\mathbf{\xi \in [1,N]}$", loc="right")
        ax.ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
        ax.title(str(f"N = {N}", r"dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = 0.3$"))

        return fig, ax

    def plotRatesLat(self, lat : Lattice, ham : Hamiltonian) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description

        Same as plotRates(), but automatically extracts information from Lattice and Hamiltonian objects.
        """

        #Extract information:
        N = lat.getN()
        d = lat.getd()
        rates = ham.getDecayRates()

        fig, ax = self.plotRates(N, d, rates)

        return fig, ax