"""
Author: tzs820, s240670, Emil Henningsen

Plotting module
"""

import matplotlib.pyplot as plt
from cycler import cycler

class Plots:

    def __init__(self):
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica",
            "axes.labelsize": 16,
            "axes.labelweight": "bold",
            "axes.titlesize": 20,
            "axes.titlelocation": "left",
            "axes.prop_cycle": cycler(color=["black", "red", "blue", "purple"]),
            "axes.grid": plt.grid(visible=True, axis="both", alpha=0.5, color="grey")
        })

    
        
