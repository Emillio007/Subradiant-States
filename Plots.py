"""
Author: tzs820, s240670, Emil Henningsen

Plotting module
"""

import matplotlib.pyplot as plt
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
        
