#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:39:20 2020

@author: jkin0004
"""

# %% codecell
# importing the required modules
from scipy.stats import gennorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import rheoSANS_fitOpt_Functions as rsf


# %% codecell
# Choose data to plot
select = [1, 1]
shear = 0
plot_exp = "36"
plot_sim = "47"

# set path and return pathnames
base = '../2D_annular_sector_extraction/'

folder_exp = base + 'py_annular_exp/'
folder_sim = base + 'py_annular_sim/'
folder_leg = base + 'LegendreFits/'
fileNames_leg = ['10wt_10ps_r0p045_legendreFit.txt',
                 # '10wt_25ps_legendreFit.txt',
                 '10wt_75ps_r0p045_legendreFit.txt',
                 '10wt_200ps_r0p045_legendreFit.txt',
                 '10wt_500ps_r0p045_legendreFit.txt']

csv_exp = folder_exp + 'annular_expData.csv'
csv_sim = folder_sim + 'annular_simData.csv'

fieldnames = ['index', 'filename', 'sample', 'shear']

if select[0] == 1:
    parameters = rsf.files_list_reduce(csv_exp, fieldnames)
    fileName_exp, sample_exp, shear_exp = rsf.files_to_reduce(parameters, plot_exp)
    data_exp = []
    data_leg = []
    for names in fileName_exp:
        data_exp.append(np.loadtxt(folder_exp + names, delimiter="  ", skiprows=1))
    for files in fileNames_leg:
        data_leg.append(np.loadtxt(base + folder_leg + files, delimiter="\t", skiprows=18))

if select[1] == 1:
    parameters = rsf.files_list_reduce(csv_sim, fieldnames)
    fileName_sim, sample_sim, shear_sim = rsf.files_to_reduce(parameters, plot_sim)
    data_sim = []
    for names in fileName_sim:
        data_sim.append(np.loadtxt(folder_sim + names, delimiter="  ", skiprows=1))

# %% codecell
colors = ["#67204C",
          "#8E177D",
          "#A91CB5",
          "#A353D9",
          "#8E83DF",
          "#86A5E4",
          "#8FC2EA",
          "#9FDEF5"]
#
# colors = ["#5B3794",
#           "#9953A1",
#           "#C87AAD",
#           "#EBA8BA",
#           "#F8DCD9"]

cols2 = ['red', 'blue']
for i in range(len(colors)):
    cols2.append(colors[len(colors)-1-i])

marker = ['o', 'D', '<', 's', '>']
# plotting the data

if select[0] == 1:
    for i, sets in enumerate(data_exp):
        c = cm.magma((i)/4., 1)  # set colour from colourmap
        # c = colors[i]
        # c = cols2[i]
        min = np.min(sets[:, 1])
        diff = 10-min
        diff = 0
        # diff = 0
        r'I(q) / $\bf{cm^{-1}}$'
        half = int(np.round(len(sets[:, 0])/2))
        plt.plot(sets[0:half, 0], sets[0:half, 2]+diff, marker=marker[i], color=c[-1],
                 markeredgewidth=0.4, linewidth=0, markersize=6, alpha=1, label='Exp')
        # shear_exp[i] + ' ' + r's$\bf{^{-1}}$'
        # plt.plot(data_leg[i][:, 0], data_leg[i][:, 2]+diff, color=c,
        #          alpha=1, linewidth=2.5)
        # plt.text(-1.5, 55, r'$\bf{\dot{\gamma} = x s^{-1}}$')

if select[1] == 1:
    for i, sets in enumerate(data_sim):
        # c = cm.cool((i+1)/5., 1)  # set colour from colourmap
        # c = colors[i]
        c = cols2[i]
        min = np.min(sets[:, 1])
        diff = 25-min
        diff = 0
        r'I(q) / $\bf{cm^{-1}}$'
        half = int(np.round(len(sets[:, 0])/2))
        plt.plot(sets[0:half, 0], sets[0:half, 1]+diff, marker='o', color='blue',
                 label='Model', markeredgecolor='blue',
                 markeredgewidth=0.4, linewidth=0, markersize=6, alpha=1)
# label=shear_sim[i] + ' ' + r's$^{-1}$'

plt.legend()
# plt.title('25 wt%', fontweight='bold')
fontsize = '16'

# naming the axes
plt.xlabel('Angle / rad', fontweight='bold', fontsize=fontsize)
plt.ylabel(r'$I \:/\: \bf{cm^{-1}}$', fontweight='bold', fontsize=fontsize)

# Set plot characteristics from rc parameters
# Axes
plt.rc('axes', linewidth=1.5)
plt.rc('axes', grid=False)
plt.rc('axes', labelsize='small')
plt.rc('axes.spines', top=False)
plt.rc('axes.spines', right=False)
#plt.rc('axes', titlesize = 'large')
#plt.rc('axes', titlelocation = 'center')

plt.xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5], [-1.5, '', '', 0, '', '', 1.5])
plt.yticks([10, 15, 20, 25, 30], [10, 15, 20, 25, 30])
plt.ylim(10, 30)
# Font
plt.rc('font', family='sans-serif')
plt.rc('font', weight='bold')
plt.rc('font', size=fontsize)

# Figure
cmToInch = 0.393701
fig_width = 15.8 * cmToInch
fig_height = 11.8 * cmToInch
plt.rc('figure', figsize=[fig_width, fig_height])
plt.rc('figure', dpi='500')
#plt.rc('figure.subplot', hspace = '0.01')
#plt.rc('figure.subplot', wspace = '0.01')

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
plt.rc('legend', loc='upper right')

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
fig.savefig('annular_10.png'.format(shear), dpi=500)

# function to show the plot
plt.show()

# %% codecell

# fig, ax = plt.subplots(1, 1)
#
# beta = 2
# x = np.linspace(gennorm.ppf(0.01, beta),
#                 gennorm.ppf(0.99, beta), 100)
# ax.plot(x, gennorm.pdf(x, beta),
#         'r-', lw=5, alpha=0.6, label='gennorm pdf')

# %% codecell
# plot I_max - I_0  vs shear

# if select[0] == 1:
#     I_ave_exp = np.mean(data_exp[0][:, 1])
#
#     shear_exp_int = []
#     I_max = []
#     for i, vals in enumerate(shear_exp):
#         shear_exp_int.append(int(vals))
#         I_max.append(np.max(data_exp[i]))
#     plt.plot(shear_exp_int, I_max, 'o-')
#
# if select[1] == 1:
#     I_ave_sim = np.mean(data_sim[0][:, 1])
#
#     shear_sim_int = []
#     I_max = []
#     for i, vals in enumerate(shear_sim):
#         shear_sim_int.append(int(vals))
#         I_max.append(np.max(data_sim[i]))
#     plt.plot(shear_sim_int, I_max)
