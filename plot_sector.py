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
plot_vert = "83,81"
plot_horiz = plot_vert

# set path and return pathnames
folder_vert = '2D_annular_sector_extraction/VertPi20/'
folder_horiz = '2D_annular_sector_extraction/HorizPi20/'

csv_filename_vert = folder_vert + 'vert_expData.csv'
csv_filename_horiz = folder_horiz + 'horiz_expData.csv'

fieldnames = ['index', 'filename', 'sample', 'shear']
        
parameters_vert = rsf.files_list_reduce(csv_filename_vert, fieldnames)
fileName_vert, sample_vert, shear_vert = rsf.files_to_reduce(parameters_vert, plot_vert)
pathName_vert = []

parameters_horiz = rsf.files_list_reduce(csv_filename_horiz, fieldnames)
fileName_horiz, sample_horiz, shear_horiz = rsf.files_to_reduce(parameters_horiz, plot_horiz)
pathName_horiz = []

vert_angle = ', 0^{\circ}}$'
horiz_angle = ', 90^{\circ}}$'

labels_vert =[]
labels_horiz = []
for i in range(len(shear_vert)):
    labels_vert.append(shear_vert[i][0:-2] + vert_angle )
    labels_horiz.append(shear_horiz[i][0:-2] + horiz_angle)

for names_vert in fileName_vert:
    pathName_vert.append(folder_vert + names_vert)
    
for names_horiz in fileName_horiz:
    pathName_horiz.append(folder_horiz + names_horiz)
    

#load data
data_vert = []
data_horiz = []

for names in pathName_vert:
     data_vert.append(np.loadtxt(names, delimiter = "\t", skiprows = 1))
     
for names in pathName_horiz:
     data_horiz.append(np.loadtxt(names, delimiter = "\t", skiprows = 1))
       
n = len(data_vert)
n=12
color=iter(cm.Paired(np.linspace(0,1,n)))

color = ['indigo', 'mediumvioletred']

# plotting the data
for i, sets in enumerate(data_vert):
#    c = cm.PRGn(i/2.,1) # set colour from colourmap
    plt.loglog(sets[:,0], sets[:,1],'o', color=color[i], label=labels_vert[i], alpha=1)
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')

# plotting the data
for i, sets in enumerate(data_horiz):
#    c = cm.PRGn(i,1) # set colour from colourmap
    plt.loglog(sets[:,0], sets[:,1],'s', color='w', markeredgecolor=color[i], label=labels_horiz[i])
#    plt.plot(sets[:,0], sets[:,1], '--', color='black')

fontsize = '10'
plt.legend()

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
plt.rc('legend', loc='lower left')


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