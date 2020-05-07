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
import rheoSANSFunctions_fitOpt_lsq as rsf
  
# Choose data to plot
plot = "6,34,44,61,71,83"

# set path and return pathnames for data
folder = '2D_annular_sector_extraction/Radially_ave/'

csv_filename = folder + 'radAve_expData.csv'
fieldnames = ['index', 'filename', 'sample', 'shear']

        
parameters = rsf.files_list_reduce(csv_filename, fieldnames)
fileName, sample, shear = rsf.files_to_reduce(parameters, plot)
pathName = []
        
for names in fileName:
    pathName.append(folder + names)

#load data
data = []
for names in pathName:
     data.append(np.loadtxt(names, delimiter = "\t", skiprows = 1))
       
     

# set path and return pathnames for sim fits

folderFits = 'radAve_fits/'  

csv_filenameFits = folderFits + 'radAve_simFits.csv'
fieldnamesFits = ['index', 'filename', 'sample', 'shear']

        
parametersFits = rsf.files_list_reduce(csv_filenameFits, fieldnames)
fileNameFits, sampleFits, shearFits = rsf.files_to_reduce(parametersFits, plot)
pathNameFits = []
        
for names in fileNameFits:
    pathNameFits.append(folderFits + names)

#load data
dataFits = []
for names in pathNameFits:
     dataFits.append(np.loadtxt(names, delimiter = "  ", skiprows = 1))   




# plotting the data
for i, sets in enumerate(data):
    c = cm.plasma(i/6.,1) # set colour from colourmap
    plt.loglog(sets[:,0], sets[:,1],'o', color=c, label=sample[i])
    plt.plot(dataFits[i][:,0], dataFits[i][:,1],'-', color='white', alpha=0.5)
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')

#for i, sets in enumerate(dataFits):
#    c = cm.plasma(i/6.,0.5) # set colour from colourmap
#    plt.plot(sets[:,0], sets[:,1],'-', color='white', alpha=0.5)
##    plt.plot(sets[:,0], sets[:,1], '--', color='black')


plt.legend()
fontsize = '10'

# naming the axes 
plt.xlabel(r'q / $\bf{\AA^{-1}}$',fontweight = 'bold', fontsize = fontsize)
plt.ylabel(r'I(q) / $\bf{cm^{-1}}$', fontweight = 'bold', fontsize = fontsize)
    
#Set plot characteristics from rc parameters
#Axes
plt.rc('axes',linewidth=1.5)
plt.rc('axes',grid=False)
plt.rc('axes', labelsize = 'small')
#plt.rc('axes', titlesize = 'large')
#plt.rc('axes', titlelocation = 'center')

#Font
plt.rc('font',family = 'sans-serif')
plt.rc('font', weight= 'bold')
plt.rc('font', size= fontsize)

#Figure
cmToInch = 0.393701
fig_width = 10.41 *cmToInch
fig_height = 13.54 *cmToInch
plt.rc('figure', figsize= [fig_width,fig_height])
plt.rc('figure', dpi= '150')
#plt.rc('figure.subplot', hspace = '0.01')
#plt.rc('figure.subplot', wspace = '0.01')

plt.rc('figure.constrained_layout', use=True)

#Grid
plt.rc('grid', color = 'b')  ## grid color
plt.rc('grid', linestyle = '-')       ## solid
plt.rc('grid', linewidth = 1) 
plt.rc('grid', alpha = 0.5) 
#grid.linewidth : 0.8     ## in points
#grid.alpha     : 1.0     ## transparency, between 0.0 and 1.0

#Legend
plt.rc('legend', frameon = False)

#Ticks
plt.rc('xtick',bottom=True)
plt.rc('ytick',left=True)
plt.rc('xtick.major',width=1.5)
plt.rc('ytick.major',width=1.5)
plt.rc('xtick.minor',width=1.5)
plt.rc('ytick.minor',width=1.5)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

  
# function to show the plot 
plt.show() 