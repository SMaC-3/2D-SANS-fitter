#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:39:20 2020

@author: jkin0004
"""

# importing the required modules
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import numpy as np
import rheoSANS_fitOpt_Functions as rsf

# Choose data to plot
select = [1, 1, 0, 0]

plot_vert_exp = "4,7"
plot_horiz_exp = "4,7"
plot_vert_sim = "4,7"
plot_horiz_sim = "4,7"

# set path and return pathnames
base = '../2D_annular_sector_extraction/'

folder_vert_exp = base + 'py_sect_vert_exp/'
folder_horiz_exp = base + 'py_sect_horiz_exp/'
folder_vert_sim = base + 'py_sect_vert_sim/'
folder_horiz_sim = base + 'py_sect_horiz_sim/'

csv_vert_exp = base + 'py_sect_vert_exp/vert_expData.csv'
csv_horiz_exp = base + 'py_sect_horiz_exp/horiz_expData.csv'
csv_vert_sim = base + 'py_sect_vert_sim/vert_simData.csv'
csv_horiz_sim = base + 'py_sect_horiz_sim/horiz_simData.csv'

fieldnames = ['index', 'filename', 'sample', 'shear']

vert_angle = ' s' + r'$^{-1}, 0^{\circ}}$'
horiz_angle = ' s' + r'$^{-1}, 90^{\circ}}$'

if select[0] == 1:
    parameters_vert_exp = rsf.files_list_reduce(csv_vert_exp, fieldnames)
    fileName_vert_exp, sample_vert_exp, shear_vert_exp = rsf.files_to_reduce(
        parameters_vert_exp, plot_vert_exp)
    data_vert_exp = []
    labels_vert_exp = []
    for i, names_vert in enumerate(fileName_vert_exp):
        labels_vert_exp.append(shear_vert_exp[i] + vert_angle)
        data_vert_exp.append(np.loadtxt(folder_vert_exp + names_vert,
                                        delimiter="  ", skiprows=1))

if select[1] == 1:
    parameters_horiz_exp = rsf.files_list_reduce(csv_horiz_exp, fieldnames)
    fileName_horiz_exp, sample_horiz_exp, shear_horiz_exp = rsf.files_to_reduce(
        parameters_horiz_exp, plot_horiz_exp)
    data_horiz_exp = []
    labels_horiz_exp = []
    for i, names_horiz in enumerate(fileName_horiz_exp):
        labels_horiz_exp.append(shear_horiz_exp[i] + horiz_angle)
        data_horiz_exp.append(np.loadtxt(folder_horiz_exp + names_horiz,
                                         delimiter="  ", skiprows=1))

if select[2] == 1:
    parameters_vert_sim = rsf.files_list_reduce(csv_vert_sim, fieldnames)
    fileName_vert_sim, sample_vert_sim, shear_vert_sim = rsf.files_to_reduce(
        parameters_vert_sim, plot_vert_sim)
    data_vert_sim = []
    labels_vert_sim = []
    for i, names_vert in enumerate(fileName_vert_sim):
        labels_vert_sim.append(shear_vert_sim[i] + vert_angle)
        data_vert_sim.append(np.loadtxt(folder_vert_sim + names_vert,
                                        delimiter="  ", skiprows=1))

if select[3] == 1:
    parameters_horiz_sim = rsf.files_list_reduce(csv_horiz_sim, fieldnames)
    fileName_horiz_sim, sample_horiz_sim, shear_horiz_sim = rsf.files_to_reduce(
        parameters_horiz_sim, plot_horiz_sim)
    data_horiz_sim = []
    labels_horiz_sim = []
    for i, names_horiz in enumerate(fileName_horiz_sim):
        labels_horiz_sim.append(shear_horiz_sim[i] + horiz_angle)
        data_horiz_sim.append(np.loadtxt(folder_horiz_sim + names_horiz,
                                         delimiter="  ", skiprows=1))


color = ['indigo', 'mediumvioletred']

# plotting the data
if select[0] == 1:
    for i, sets in enumerate(data_vert_exp):
        plt.loglog(sets[:, 0], sets[:, 1], 'o', color=color[i],
                   label=labels_vert_exp[i], alpha=1)
if select[1] == 1:
    for i, sets in enumerate(data_horiz_exp):
        plt.loglog(sets[:, 0], sets[:, 1], 's', color='w',
                   markeredgecolor=color[i], label=labels_horiz_exp[i])

if select[2] == 1:
    for i, sets in enumerate(data_vert_sim):
        plt.loglog(sets[:, 0], sets[:, 1], 'o', color=color[i],
                   label=labels_vert_sim[i], alpha=1)

if select[3] == 1:
    for i, sets in enumerate(data_horiz_sim):
        plt.loglog(sets[:, 0], sets[:, 1], 's', color='w',
                   markeredgecolor=color[i], label=labels_horiz_sim[i])


# Modify plot properties
fontsize = '10'
plt.legend()

# naming the axes
plt.xlabel(r'q / $\bf{\AA^{-1}}$', fontweight='bold', fontsize=fontsize)
plt.ylabel(r'I(q) / $\bf{cm^{-1}}$', fontweight='bold', fontsize=fontsize)

# Set plot characteristics from rc parameters
# Axes
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
plt.rc('axes', linewidth=1.5)
plt.rc('axes', grid=False)
plt.rc('axes', labelsize='small')
#plt.rc('axes', titlesize = 'large')
#plt.rc('axes', titlelocation = 'center')

# plt.xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5], [-1.5, '', '', 0, '', '', 1.5])
# plt.yticks([30, 40, 50, 60, 70], [30, '', 50, '', 70])

# Font
plt.rc('font', family='sans-serif')
plt.rc('font', weight='bold')
plt.rc('font', size=fontsize)

# Figure
cmToInch = 0.393701
fig_width = 10.41 * cmToInch
fig_height = 13.54 * cmToInch
plt.rc('figure', figsize=[fig_width, fig_height])
plt.rc('figure', dpi='150')
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
plt.rc('legend', loc='lower left')

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
fig.savefig('sectors25wt.png', dpi=300)

# function to show the plot
plt.show()
