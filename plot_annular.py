#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:39:20 2020

@author: jkin0004
"""

# importing the required modules
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import rheoSANS_fitOpt_Functions as rsf

# Choose data to plot
plot = "13,15,16,18,21"

# set path and return pathnames
folder = '../2D_annular_sector_extraction/py_annular/'

csv_filename = folder + 'annular_expData.csv'
fieldnames = ['index', 'filename', 'sample', 'shear']

parameters = rsf.files_list_reduce(csv_filename, fieldnames)
fileName, sample, shear = rsf.files_to_reduce(parameters, plot)
pathName = []
print(shear)

# fileName = np.insert(fileName, 0, fileName[-1])
# sample = np.insert(sample, 0, sample[-1])
# shear = np.insert(shear, 0, shear[-1])
#
# fileName = fileName[0:-1]
# sample = sample[0:-1]
# shear = shear[0:-1]

# labels = 'iii'

for names in fileName:
    pathName.append(folder + names)

# load data
data = []
for names in pathName:
    data.append(np.loadtxt(names, delimiter="  ", skiprows=1))

colors = ["#67204C",
          "#8E177D",
          "#A91CB5",
          "#A353D9",
          "#8E83DF",
          "#86A5E4",
          "#8FC2EA",
          "#9FDEF5"]

colors = ["#5B3794",
          "#9953A1",
          "#C87AAD",
          "#EBA8BA",
          "#F8DCD9"]

cols2 = []
for i in range(len(colors)):
    cols2.append(colors[len(colors)-1-i])

print(cols2)

# marker = ['o', 'D', '<', '*', '>']
# plotting the data
for i, sets in enumerate(data):
    # c = cm.cool((i+1)/5., 1)  # set colour from colourmap
    # c = colors[i]
    c = cols2[i]
    min = np.min(sets[:, 1])
    diff = 25-min
    r'I(q) / $\bf{cm^{-1}}$'
    plt.plot(sets[0:100, 0], sets[0:100, 1]+diff, marker='o', color=c, label=shear[i] + ' ' + r's$^{-1}$', markeredgecolor='black',
             markeredgewidth=0.4, linewidth=0, markersize=6)
    # plt.plot(sets[:, 0], sets[:, 1]+diff, 'o', color=c, label=labels)
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')


plt.legend()
# plt.title('25 wt%', fontweight='bold')
fontsize = '10'

# naming the axes
plt.xlabel('Angle / rad', fontweight='bold', fontsize=fontsize)
plt.ylabel(r'$I \:/\: \bf{cm^{-1}}$', fontweight='bold', fontsize=fontsize)

# Set plot characteristics from rc parameters
# Axes
# plt.rc('axes', linewidth=1.5)
# plt.rc('axes', grid=False)
# plt.rc('axes', labelsize='small')
# # plt.rc('axes', titlesize = 'large')
# # plt.rc('axes', titlelocation = 'center')
#
# # Font
# plt.rc('font', family='sans-serif')
# plt.rc('font', weight='bold')
# plt.rc('font', size=fontsize)
#
# # Figure
# cmToInch = 0.393701
# fig_width = 15.41 * cmToInch
# fig_height = 13.54 * cmToInch
# plt.rc('figure', figsize=[fig_width, fig_height])
# plt.rc('figure', dpi='150')
# # plt.rc('figure.subplot', hspace = '0.01')
# # plt.rc('figure.subplot', wspace = '0.01')
#
# plt.rc('figure.constrained_layout', use=True)
#
# # Grid
# plt.rc('grid', color='b')  # grid color
# plt.rc('grid', linestyle='-')  # solid
# plt.rc('grid', linewidth=1)
# plt.rc('grid', alpha=0.5)
# # grid.linewidth : 0.8     ## in points
# # grid.alpha     : 1.0     ## transparency, between 0.0 and 1.0
#
# # Legend
# plt.rc('legend', frameon=False)
#
# # Ticks
# plt.rc('xtick', bottom=True)
# plt.rc('ytick', left=True)
# plt.rc('xtick.major', width=1.5)
# plt.rc('ytick.major', width=1.5)
# plt.rc('xtick.minor', width=1.5)
# plt.rc('ytick.minor', width=1.5)
# plt.xticks(fontsize=fontsize)
# plt.yticks(fontsize=fontsize)
#
#
# # function to show the plot
# plt.show()

# Set plot characteristics from rc parameters
# Axes
plt.rc('axes', linewidth=1.5)
plt.rc('axes', grid=False)
plt.rc('axes', labelsize='small')
#plt.rc('axes', titlesize = 'large')
#plt.rc('axes', titlelocation = 'center')

plt.xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5], [-1.5, '', '', 0, '', '', 1.5])
plt.yticks([30, 40, 50, 60, 70], [30, '', 50, '', 70])

# Font
plt.rc('font', family='sans-serif')
plt.rc('font', weight='bold')
plt.rc('font', size=fontsize)

# Figure
cmToInch = 0.393701
fig_width = 13.54 * cmToInch
fig_height = 10.41 * cmToInch
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
fig.savefig('annular25exp.png', dpi=300)

# function to show the plot
plt.show()
