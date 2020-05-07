#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:08:51 2020

@author: jkin0004
"""

import rheoSANSFunctions_fitOpt_lsq as rsf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

indexSelect = '94'

indexNums = rsf.evaluate_files_list(indexSelect)

sans = rsf.sans2d()
sans.qmin = 0.07
grids = []

q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
xxq, yyq = np.meshgrid(q,q)




for idx in indexNums:
    sans.getData(str(idx))
    interp, zzq = sans.interpData(sans.expData)
    grids.append(zzq)
#    sans.getSim(str(idx))
#    interp, zzq = sans.interpData(sans.simImport)
#    grids.append(zzq)
    
    del interp, zzq, sans.expData

minima = []
maxima = []

for sets in grids:
    minima.append(np.nanmin(sets))
    maxima.append(np.nanmax(sets))

vmin = np.amin(minima)
vmax = np.amax(maxima)

cmToInch = 0.393701

for i in range(len(indexNums)):
    plt.figure(figsize = [5.1, 4.08], dpi = 200)
    plt.pcolormesh(xxq, yyq, grids[i], cmap='jet', vmin=vmin, vmax=vmax, norm=matplotlib.colors.LogNorm())
    plt.title(indexNums[i])
    # set the limits of the plot to the limits of the data
    #plt.axis([sim.xxq.min(), sim.xxq.max(), sim.yyq.min(), sim.yyq.max()])
    plt.colorbar()
    
    fontsize = '10'

    # naming the axes 
    plt.xlabel(r'$q_x$',fontweight = 'normal', fontsize = fontsize)
    plt.ylabel(r'$q_y$', fontweight = 'normal', fontsize = fontsize)
        
    #Set plot characteristics from rc parameters
    #Axes
    plt.rc('axes',linewidth=1)
    plt.rc('axes',grid=False)
    plt.rc('axes', labelsize = 'small')
    #plt.rc('axes', titlesize = 'large')
    #plt.rc('axes', titlelocation = 'center')
    
    #Font
    plt.rc('font',family = 'sans-serif')
    plt.rc('font', weight= 'normal')
    plt.rc('font', size= fontsize)
    
    #Figure
    
    #fig_width = 10 *cmToInch
    #fig_height = 10 *cmToInch
    #plt.rc('figure', figsize= [fig_width,fig_height])
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
#    plt.xticks(rotation='vertical')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    
    


    
    

    