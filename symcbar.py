"""Credit: Nicolas, TA in Continuum Mechanics 2023, UCPH"""

import matplotlib as mpl
cmap = mpl.cm.RdBu; # red-blue divergent colormap
clim = 0.5 # colorbar max/min range
ticks=np.linspace(-clim,clim,10);
# ... calc tauxx here ...
h=plot(tauxx*1e-6, title='tau_xx [MPa]', cmap=cmap, norm=mpl.colors.Normalize(vmin=ticks[0], vmax=ticks[-1]))
hcb=plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=ticks[0], vmax=ticks[-1]), cmap=cmap),ticks=ticks)
