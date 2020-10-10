#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:08:51 2020

@author: jkin0004
"""

import rheoSANS_fitOpt_Functions as rsf
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.axes as ax
import matplotlib.colors as mcolors
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import nmmn.plots
wolfram = nmmn.plots.wolframcmap()  # for Mathematica's cmap
parula = nmmn.plots.parulacmap()  # for MATLAB's cmap
# turbo=nmmn.plots.turbocmap() # Turbo


expInd = '41'
simInd = '53'

indexNums = rsf.evaluate_files_list(expInd)
# indexNums2 = rsf.evaluate_files_list(simInd)

sans = rsf.sans2d()
sans.qmin = 0.007
sans.nq = 500
sans.q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
sans.xxq, sans.yyq = np.meshgrid(sans.q, sans.q)
grids = []
plot = []

q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
xxq, yyq = np.meshgrid(q, q)

for i, idx in enumerate(indexNums):
    sans.getData(str(idx))
    # sans.getSim(str(indexNums2[i]))
    interp_exp, zzq_exp = sans.interpData(data=None, dataSet=sans.expData)
    # interp_sim, zzq_sim = sans.interpData(data=None, dataSet=sans.simImport)
    # zzq = abs(zzq_exp - zzq_sim)
    grids.append(zzq_exp)
    plot.append('exp' + str(idx))
    # grids.append(zzq_sim)
    # plot.append('sim' + str(indexNums2[i]))
    # grids.append(zzq)

    # sans.getSim(str(idx))
    # interp, zzq = sans.interpData(data=None, dataSet=sans.simImport)
    # grids.append(zzq)

    # del interp, zzq_exp, sans.expData

minima = []
maxima = []

for sets in grids:
    minima.append(np.nanmin(sets[sets > 0]))
    maxima.append(np.nanmax(sets))

vmin = np.amin(minima)
vmax = np.amax(maxima)
print(vmax)

cmToInch = 0.393701

# colors1 = plt.cm.Blues(np.linspace(0.3, 1, 500))
# colors2 = plt.cm.jet(np.linspace(0, 1, 500))
#
# # combine them and build a new colormap
# colors = np.vstack((colors1, colors2))
# mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
# print(vmin)

# plot annular ring
circ = np.linspace(0, 2*np.pi, 50)
rad = 0.045
thx = 0.01
inner_x = rad*np.cos(circ)
inner_y = rad*np.sin(circ)
outer_x = (thx+rad)*np.sin(circ)
outer_y = (thx+rad)*np.cos(circ)
cmap = parula

for i in range(len(plot)):

    plt.pcolormesh(xxq, yyq, grids[i], cmap=cmap, vmin=2,
                   vmax=vmax, norm=mcolors.LogNorm())

    plt.plot(inner_x, inner_y, 'black', lineStyle='-')
    plt.plot(outer_x, outer_y, 'black', lineStyle='-')

    # plt.plot([-0.045, 0.045], [-0.19, 0.19], color='black')
    # plt.plot([-0.045, 0.045], [0.19, -0.19], color='black')
    plt.title(plot[i])
    # set the limits of the plot to the limits of the data
    # plt.axis([sim.xxq.min(), sim.xxq.max(), sim.yyq.min(), sim.yyq.max()])
    plt.colorbar()

    fontsize = '10'

    # naming the axes
    plt.xlabel(r'$q_x$', fontweight='normal', fontsize=fontsize)
    plt.ylabel(r'$q_y$', fontweight='normal', fontsize=fontsize)

    # Set plot characteristics from rc parameters
    # Axes
    plt.rc('axes', linewidth=1)
    plt.rc('axes', grid=False)
    plt.rc('axes', labelsize='small')
    ax = plt.gca()
    ax.set_aspect(aspect='equal')
    # plt.rc('axes', titlesize = 'large')
    # plt.rc('axes', titlelocation = 'center')

    # Font
    plt.rc('font', family='sans-serif')
    plt.rc('font', weight='bold')
    plt.rc('font', size=fontsize)

    # Figure

    # fig_width = 10 *cmToInch
    # fig_height = 10 *cmToInch
    # plt.rc('figure', figsize= [fig_width,fig_height])
    plt.rc('figure', dpi='500')
    # plt.rc('figure.subplot', hspace = '0.01')
    # plt.rc('figure.subplot', wspace = '0.01')

    plt.rc('figure.constrained_layout', use=True)

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
    plt.rc('ytick', left=False)
    plt.xticks(rotation=20)
    # plt.yticks(rotation=45)
    plt.rc('xtick.major', width=1.5)
    plt.rc('ytick.major', width=1.5)
    # plt.rc('xtick.minor', width=1.5)
    # plt.rc('ytick.minor', width=1.5)
#    plt.xticks(rotation='vertical')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig('10wt%_500ps.png'.format(i), dpi=500)


# %% codecell
# plt.rcParams["figure.figsize"] = 12.8, 9.6
# # Normalize the colors based on Z value
# zzq = grids[0]
# norm = plt.Normalize(zzq.min(), zzq.max())
# colors = cm.jet(norm(zzq))
# ax = plt.axes(projection='3d')
# zzq1 = np.ones_like(zzq)
# neg = xxq < 0
# zzq1 = zzq[~neg]
# # zzq[0:100] = 1
# surf = ax.plot_surface(xxq, yyq, zzq, cmap=cm.jet,
#                        linewidth=0, antialiased=False,
#                        shade=False, rstride=1, cstride=1)
# surf.set_facecolor((0, 0, 0, 0))
