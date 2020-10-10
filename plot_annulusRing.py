# %% codecell
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import sasmodels.compare
import sasmodels.core
import sasmodels.direct_model
import sasmodels.data
import rheoSANS_fitOpt_Functions as rsf
import annular_sector_extraction_loop as ansect
import itertools

import nmmn.plots
wolfram = nmmn.plots.wolframcmap()  # for Mathematica's cmap
parula = nmmn.plots.parulacmap()  # for MATLAB's cmap
# turbo=nmmn.plots.turbocmap() # Turbo

# %% codecell
sans = rsf.sans2d()
sans.nq = 500
sans.qmin = 0.007
sans.q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
sans.xxq, sans.yyq = np.meshgrid(sans.q, sans.q)

source = 'exp'

if source == 'sim':
    sans.getSim('45')
    interp, a = sans.interpData_noMask(None, sans.simImport)
elif source == 'exp':
    sans.getData('41')
    interp, a = sans.interpData_noMask(None, sans.expData)

bins, I, err, annul_data = ansect.annular(interp, 0.045, 0.01)
annul_x = annul_data[0]
annul_y = annul_data[1]
annul_I = annul_data[2]
I_ann = annul_data[3]

# %% codecell
interp.data[~I_ann] = np.nan


a = np.reshape(interp.qx_data, (sans.nq, sans.nq))
b = np.reshape(interp.qy_data, (sans.nq, sans.nq))
c = np.reshape(interp.data, (sans.nq, sans.nq))
# vmax = np.nanmax(c)
# vmin = np.nanmin(c)
vmax = np.nanmax(sans.expData.data)
vmin = np.nanmin(sans.expData.data)
print(vmax)

# d = np.reshape(I_ann, (50,50))
# ax = plt.subplot()
plt.figure(figsize=[5, 5], dpi=500)
# plt.plot(np.array([0, 0.06]), [0, 0.06])
plt.pcolormesh(a, b, c, cmap=parula, vmin=2, vmax=vmax, norm=mcolors.LogNorm())
# plt.colorbar(shrink=0.7, ticks=[20, 25, 30, 35])
plt.xlim(-0.055, 0.055)
plt.ylim(-0.055, 0.055)
plt.xticks([-0.05, -0.025, 0, 0.025, 0.05], [-50, -25, 0, 25, 50])
plt.yticks([-0.05, -0.025, 0, 0.025, 0.05], [-50, -25, 0, 25, 50])
fontsize = '16'

plt.xlabel('q ' + r'/ $\bf{10^{-3} \: \: \AA^{-1}}$', fontweight='bold', fontsize=fontsize)
plt.ylabel('q ' + r'/ $\bf{10^{-3} \: \: \AA^{-1}}$', fontweight='bold', fontsize=fontsize)
plt.rc('figure', frameon=False)
plt.rc('axes.spines', top=False)
plt.rc('axes.spines', right=False)
plt.rc('axes.spines', left=True)
plt.rc('axes.spines', bottom=True)
plt.rc('ytick', left=True)

# Set plot characteristics from rc parameters
# Axes
ax = plt.gca()
ax.set_aspect(aspect='equal')
plt.rc('axes', linewidth=1.5)
plt.rc('axes', grid=False)
plt.rc('axes', labelsize='small')
# plt.rc('axes', titlesize = 'large')
# plt.rc('axes', titlelocation = 'center')

# Font
plt.rc('font', family='sans-serif')
plt.rc('font', weight='bold')
plt.rc('font', size=fontsize)

# Figure
# cmToInch = 0.393701
# fig_width = 10.41 * cmToInch
# fig_height = 13.54 * cmToInch
# plt.rc('figure', figsize=[fig_width, fig_height])
# plt.rc('figure', dpi='150')
# plt.rc('figure.subplot', hspace = '0.01')
# plt.rc('figure.subplot', wspace = '0.01')

# plt.rc('figure.constrained_layout', use=True)

# Grid
plt.rc('grid', color='b')  # grid color
plt.rc('grid', linestyle='-')  # solid
plt.rc('grid', linewidth=1)
plt.rc('grid', alpha=0.5)
# grid.linewidth : 0.8     ## in points
# grid.alpha     : 1.0     ## transparency, between 0.0 and 1.0

# Legend
plt.rc('legend', frameon=False)

# Ticks
plt.rc('xtick', bottom=True)
plt.rc('ytick', left=True)
plt.rc('xtick.major', width=1.5)
plt.rc('ytick.major', width=1.5)
plt.rc('xtick.minor', width=1.5)
plt.rc('ytick.minor', width=1.5)
plt.rc('xtick.major', size=6)
plt.rc('ytick.major', size=6)
plt.rc('xtick.minor', size=4)
plt.rc('ytick.minor', size=4)

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

fig = plt.figure(1)
fig.savefig('annularRing_10wt_500ps.png', dpi=500)

# %% codecell

# plt.rcParams["figure.figsize"] = 12.8, 9.6
# # Normalize the colors based on Z value
# norm = plt.Normalize(zzq.min(), zzq.max())
# colors = cm.jet(norm(zzq))
# ax = plt.axes(projection='3d')
# surf = ax.plot_surface(self.xxq, self.yyq, zzq, facecolors=colors, shade=False)
# surf.set_facecolor((0, 0, 0, 0))
