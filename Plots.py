"""
Author: tzs820, s240670, Emil Henningsen

Plotting module
"""

import matplotlib.pyplot as plt
import matplotlib
import Lattice
import Hamiltonian
from numpy import ndarray
from cycler import cycler
from typing import Literal

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

    def plotDipolesPlane(self, lat : Lattice, ax : plt.Axes = None, plane : Literal["xy", "xz", "yz"] = "xy", title : str = None, 
                         xlim : tuple[float, float] = None, ylim : tuple[float, float] = None, ham : Hamiltonian = None, 
                         index : int = None, legend : bool = True) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description
        """
        from numpy import array
        firstTime = False
        if ax is None:
            plt.figure()
            ax = plt.gca()
            firstTime = True    #Raise flag for colorbar
        fig = plt.gcf()

        ampl = array([1])  #just the same color, if no ham provided
        if not ham is None and not index is None:
            ampl = ham.getAmplNorm(index=index)

        pos, dir, pola = lat.getPositions(), lat.getDisplacements(), lat.getPolarizations()

        x = pos[:,0]    #All x coordinates
        y = pos[:,1]    #All y coordinates
        z = pos[:,2]    #All z coordinates

        polax = pola[:,0]
        polay = pola[:,1]
        polaz = pola[:,2]

        ax1, ax2 = None, None
        p1, p2 = None, None
        po1, po2 = None, None
        match(plane):
            case "xy":
                ax1, ax2 = "x", "y"
                p1, p2 = x, y
                po1, po2 = polax, polay
            case "xz":
                ax1, ax2 = "x", "z"
                p1, p2 = x, z
                po1, po2 = polax, polaz
            case "yz":
                ax1, ax2 = "y", "z"
                p1, p2 = y, z
                po1, po2 = polay, polaz

        cmap = matplotlib.cm.coolwarm
        norm = matplotlib.colors.Normalize(vmin=min(ampl), vmax=max(ampl))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        ax.scatter(p1, p2, c=ampl, cmap=cmap, norm=norm, label="sites")
        ax.quiver(p1, p2, po1, po2, scale=15, width=0.005, color=cmap(norm(ampl)), label=r"$\hat{d}$", pivot="mid")
        if firstTime:
            plt.colorbar(sm, label=r"Amplitude norm of $|e_j> = |c_j|$")
        if not ylim is None:
            lower, upper = ylim
            ax.set_ylim(lower, upper)
        if not xlim is None:
            lower, upper = xlim
            ax.set_xlim(lower, upper)
        ax.set_xlabel(r"$\mathbf{\hat{%c}}$" % ax1, loc="right")
        ax.set_ylabel(r"$\mathbf{\hat{%c}}$" % ax2, loc="top")
        if title is None:
            title = r"Linear lattice of $N=50$ dipoles, $\frac{d}{\lambda_0}=0.3$"
        fig.suptitle(title, wrap = True, size="large")
        if legend:
            plt.legend()
        
        return fig, ax
    
    #Plot rates manually
    def plotRates(self, N : int, d : float, rates : ndarray, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None) -> plt.Figure:
        """
        TODO: Description
        """
        from textwrap import wrap

        if title == None: 
            title = "\n".join(wrap(r"$N = %s$ dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = %s$" % (N, d), 60))

        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        else:
            fig = ax.get_figure()
        
        ax.plot(range(1, N+1), rates, 'o')
        ax.set_yscale(scaley)
        ax.set_xscale(scalex)
        ax.set_xlabel(r"$\mathbf{\xi \in [1,%s]}$" % N, loc="right")
        ax.set_ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
        ax.set_title(title, loc="center")

        return fig

    def plotRatesLat(self, lat : Lattice, ham : Hamiltonian, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None) -> plt.Figure:
        """
        TODO: Description

        Same as plotRates(), but automatically extracts information from Lattice and Hamiltonian objects.
        """

        #Extract information:
        N = lat.getN()
        d = lat.getd()
        rates = ham.getDecayRates()

        fig = self.plotRates(N, d, rates, ax, scalex, scaley, title)

        return fig