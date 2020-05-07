#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:38:36 2020

@author: jkin0004
"""

import sans1d
import numpy as np

from scipy.optimize import least_squares

#model.parameters.kernel_parameters

#cylinder
pars2sim = (dict(scale=1,
            background=0.38,
            sld= -0.4,
            sld_solvent=6.3, 
            radius=20.3,
            radius_pd=.1, 
            radius_pd_n=35,
#            radius_cap = 20,
            length=96.8,
            length_pd=0,
            length_pd_n = 35,
            radius_effective_mode=0,
            radius_effective=36,
            volfraction=0.32,
            charge=26.9, 
            temperature=298.0,
            concentration_salt=0.38, 
            dielectconst=80.2))


indexSelected = '9'
p_list = ['background', 'length','scale', 'radius_effective','charge']
p_guess = [0.45, 73.3,0.91,32.7,26.6]
#p_bounds = ([0,200],[0,200],[1000,2000])
p_bounds = []
#p_bounds = ([25,0.25,90,18.5,0.9],[50,0.45,110,20.5,1.2])

options = [pars2sim, p_list, p_guess, p_bounds, indexSelected] 

def sansFit(options):

    pars2sim = options[0] #Update default pars with these
    p_list = options[1]
    p_guess = options[2]
    p_bounds = options[3]
    indexSelect = options[4]
    
    sans = sans1d.Sans1d()
    sans.qmin = 0.04
    sans.getData(indexSelect)
    sans.pars.update(pars2sim)
    sans.buildSimData()
    sans.makeCalc(sans.simData)
    
    optim = least_squares(sans.objFunc, p_guess, ftol = 1e-15, xtol = 1e-15, method = 'lm', args=(p_list,))
    
#    optim = least_squares(sans.objFunc, p_guess, method = 'lm', 
#                          args=(p_list,))
#    
    print(optim.x)
    out = optim
     
    sans.pars.update(dict(zip(p_list,optim.x)))
    optimSim = sans.calculator(**sans.pars)
    
    sans.sasPlot(data = sans.expData, sim = optimSim)  
        
    fileType = '.dat'
    fileType2 = '.txt'
    
    save = input('Would you like to save the data? Input 1 for yes [y]: ')
    
    if save == '1' or save == '':
        location = '../1D_simFits/ReducedChi2_fits/'
        sans1d.write_3_column(location + sans.description + '_sim' + fileType, sans.expData.x, optimSim)
        sans1d.saveStats(location + sans.description + '_stats' + fileType2, minParams=sans.pars, minPars=dict(zip(p_list,optim.x)), chi2='chi2_r: '+ str(np.nansum(optim.fun)), res='')
    
    return out

out = sansFit(options)