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
plot = "71,63,64,65,66,67"

# set path and return pathnames
folder = '../2D_annular_sector_extraction/Annular/'

csv_filename = folder + 'annular_expData.csv'
fieldnames = ['index', 'filename', 'sample', 'shear']
        
parameters = rsf.files_list_reduce(csv_filename, fieldnames)
fileName, sample, shear = rsf.files_to_reduce(parameters, plot)
pathName = []

fileName = np.insert(fileName, 0, fileName[-1])
sample = np.insert(sample, 0, sample[-1]) 
shear = np.insert(shear, 0, shear[-1])

fileName = fileName[0:-1]
sample = sample[0:-1]
shear = shear[0:-1]

labels = shear
        
for names in fileName:
    pathName.append(folder + names)

#load data
data = []
for names in pathName:
     data.append(np.loadtxt(names, delimiter = "\t", skiprows = 1))
     
    
       
# plotting the data
for i, sets in enumerate(data):
    c = cm.plasma(i/6.,1) # set colour from colourmap
    plt.plot(sets[:,0], sets[:,1],'o', color=c, label=labels[i])
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')


plt.legend()
plt.title('20 wt%', fontweight = 'bold')
fontsize = '10'

# naming the axes 
plt.xlabel('Angle / rad',fontweight = 'bold', fontsize = fontsize)
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