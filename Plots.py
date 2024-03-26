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

    def add(self, fig, name : str = None) -> None:
        """
        Add a figure to the instantiated manager container with given name (dict key). 
        """
        if name == None:
            number = len(self.figures)
            name = number
        self.figures[name] = fig

    def remove(self, name : str) -> None:
        """
        Remove specified figure from container (name = dict key)
        """
        self.figures.pop(name)

    def show(self) -> None:
        plt.show()

    """SET modules:"""

    def setInteractive(self, interactive : bool = False) -> None:
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

    def plot(self, *args, **kwargs) -> list:
        return plt.plot(*args, **kwargs)

    def plotDipoles(self, lat : Lattice) -> plt.Figure:
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

        fig = plt.figure()
        plt.plot(x, z, 'o', color="black", label="sites")
        plt.quiver(x, z, polax, polaz, scale=15, width=0.005, color="red", label=r"$\hat{d}$", pivot="mid")
        plt.ylim(-1, 1)
        plt.xlabel(r"$\mathbf{\hat{x}}$", loc="right")
        plt.ylabel(r"$\mathbf{\hat{z}}$", loc="top")
        plt.title(r"Linear lattice of $N=50$ dipoles, $\frac{d}{\lambda_0}=0.3$", wrap = True)
        plt.legend()
        
        return fig
    
    #Plot rates manually
    def plotRates(self, N : int, d : float, rates : ndarray, scalex : str = "linear", scaley : str = "linear", title : str = None) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description
        """
        from textwrap import wrap

        if title == None: 
            title = "\n".join(wrap(r"$N = $" + "{}".format(N) + r" dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = $" + f"{d}", 60))

        fig = plt.figure()
        plt.plot(range(1, N+1), rates, 'o')
        plt.yscale(scaley)
        plt.xscale(scalex)
        plt.xlabel(r"$\mathbf{\xi \in [1,}$" + f"{N}]", loc="right")
        plt.ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
        plt.title(title)

        return fig

    def plotRatesLat(self, lat : Lattice, ham : Hamiltonian, scalex : str = "linear", scaley : str = "linear", title : str = None) -> plt.Figure:
        """
        TODO: Description

        Same as plotRates(), but automatically extracts information from Lattice and Hamiltonian objects.
        """

        #Extract information:
        N = lat.getN()
        d = lat.getd()
        rates = ham.getDecayRates()

        fig = self.plotRates(N, d, rates, scalex, scaley, title)

        return fig