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
plot = "2, 4-8"

# set path and return pathnames for data
folder = '../2D_annular_sector_extraction/py_sect_radAve/'

csv_filename = folder + 'radAve_expData.csv'
fieldnames = ['index', 'filename', 'sample', 'shear']


parameters = rsf.files_list_reduce(csv_filename, fieldnames)
fileName, sample, shear = rsf.files_to_reduce(parameters, plot)
pathName = []

for names in fileName:
    pathName.append(folder + names)

# load data
data = []
for names in pathName:
    data.append(np.loadtxt(names, delimiter="  ", skiprows=1))


# set path and return pathnames for sim fits

folderFits = '../1D_simFits/ReducedChi2_fits/'

csv_filenameFits = folderFits + 'radAve_simFits.csv'
fieldnamesFits = ['index', 'filename', 'sample', 'shear']


parametersFits = rsf.files_list_reduce(csv_filenameFits, fieldnames)
fileNameFits, sampleFits, shearFits = rsf.files_to_reduce(parametersFits, plot)
pathNameFits = []

for names in fileNameFits:
    pathNameFits.append(folderFits + names)

# load data
dataFits = []
for names in pathNameFits:
    dataFits.append(np.loadtxt(names, delimiter="  ", skiprows=1))


# plotting the data
for i, sets in enumerate(data):
    c = cm.viridis_r((i+1)/6., 1)  # set colour from colourmap
    plt.loglog(sets[:, 0], sets[:, 1]*3**i, 'o', color=c, label=sample[i])
    plt.plot(dataFits[i][:, 0], dataFits[i][:, 1]*3**i, '-', color='red', alpha=1)
    plt.plot([0.04, 0.04], [1, 10000], 'black', lineStyle='--')
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')

# for i, sets in enumerate(dataFits):
#    c = cm.plasma(i/6.,0.5) # set colour from colourmap
#    plt.plot(sets[:,0], sets[:,1],'-', color='white', alpha=0.5)
##    plt.plot(sets[:,0], sets[:,1], '--', color='black')


plt.legend()
fontsize = '10'
plt.rc('text', usetex=False)
# plt.rc('font', family='serif')
# naming the axes
text = 'q ' + r'/ $\bf{\AA^{-1}}$'
plt.xlabel(text, fontweight='bold', fontsize=fontsize)
plt.ylabel(r'I(q) / $\bf{cm^{-1}}$', fontweight='bold', fontsize=fontsize)

ax = plt.gca()

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
# a, b = plt.yticks()
# plt.yticks([1,10,100,1000])
# ax.set_ylabels('linear')
# ax.set_yscale('linear')
# print(a)
# print(list(b))

# Set plot characteristics from rc parameters
# Axes
plt.rc('axes', linewidth=1.5)
plt.rc('axes', grid=False)
plt.rc('axes', labelsize='small')
#plt.rc('axes', titlesize = 'large')
#plt.rc('axes', titlelocation = 'center')

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
fig.savefig('SANS1D.png', dpi=300)

# function to show the plot
plt.show()
# plt.savefig('testSave.png', dpi=150)
