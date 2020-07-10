#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 16:39:54 2020

@author: jkin0004
"""

import sans1d


pars2sim = (dict(scale=1,
            background=0.23,
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
            radius_effective=32.6,
            volfraction=0.32,
            charge=26.6, 
            temperature=298.0,
            concentration_salt=0.38, 
            dielectconst=80.2))

fitChoose = dict(scale = 1,
                 sld = 0,
               background = 1,
               radius = 0,
               radius_pd = 0,
               length = 1,
               radius_effective = 0,
               charge = 1)

#fitChoose = dict(scale = 1,
#               background = 0,
#               radius = 1,
#               radius_pd = 0,
#               length = 1,
#               radius_effective = 1,
#               charge = 1)

#REMEMBER TO CHANGE THESE
indexSelected = '9'
save = 0 #Set to 1 to save

fitList = []

# fitting parameters
if fitChoose['scale'] == 1:
    par = 'scale'
    scale_min = 0.15
    scale_max = 0.18
    scale_points = 5
    
    fitList.append([[par], [scale_min, scale_max, scale_points] ])
    
if fitChoose['sld'] ==1:
    par = 'sld'
    sld_min = -0.5
    sld_max = 0
    sld_points = 5
    
    fitList.append([[par], [sld_min, sld_max, sld_points] ])

if fitChoose['background'] == 1:
    par = 'background'
    background_min = 0.2
    background_max = 0.5
    background_points = 10
    
    fitList.append([[par], [background_min, background_max, background_points] ])
    
if fitChoose['radius'] == 1:
    par = 'radius'
    radius_min = 18.5
    radius_max = 20.5
    radius_points = 10
    
    fitList.append([[par], [radius_min, radius_max, radius_points] ])
    
if fitChoose['radius_pd'] == 1:
    par = 'radius_pd'
    radius_pd_min = 0.18
    radius_pd_max = 0.19
    radius_pd_points = 1
    
    fitList.append([[par], [radius_pd_min, radius_pd_max, radius_pd_points] ])
    
if fitChoose['length'] == 1:
    par = 'length'
    length_min = 80
    length_max = 90
    length_points = 10
    
    fitList.append([[par], [ length_min, length_max, length_points] ])
        
if fitChoose['radius_effective'] == 1:
    par = 'radius_effective'
    radius_effective_min = 30
    radius_effective_max = 40
    radius_effective_points = 10
    
    fitList.append([[par], [radius_effective_min ,radius_effective_max , radius_effective_points] ])
    
if fitChoose['charge'] == 1:
    par = 'charge'
    charge_min = 25
    charge_max = 30
    charge_points = 9
    
    fitList.append([[par], [ charge_min, charge_max, charge_points] ])


#fitList = [ [ ['length'] , [80, 110, 3]   ] , [ ['radius'] , [16, 22, 4]   ], [ ['radius_pd'] , [0.1, 0.2, 1]   ] ]# , [['radius'], [19, 20, 2]] , [ ['radius_effective'], [36,39,3] ], [ ['charge'] , [30,40,3]] ] #, [['scale'] ,[1,1.1,3]] ]


options = [pars2sim, fitList, indexSelected] 

def sansFit(options, save):

    pars2sim = options[0] #Update default pars with these
    listInfo = options[1]
    indexSelect = options[2]
    
    sans = sans1d.Sans1d()
    sans.qmin = 0.0
    sans.getData(indexSelect)
    sans.pars.update(pars2sim)
    sans.buildSimData()
    sans.makeCalc(sans.simData)

    parsUpdate, paraInfo = sans1d.buildPars(listInfo)
    scatloop = sans1d.simCalc(parsUpdate, sans) #Rate determining step
    print('error min')
    listInfo, stats, minPars = sans1d.errorCalc(listInfo, scatloop, parsUpdate, paraInfo, sans)        
      
    minSim = scatloop[stats[1][1]]
    
    sans.sasPlot(data = sans.expData, sim = minSim)
    
    chi2 = 'chi2: ' + str(stats[1][0][stats[1][1]]) #chi2. 
    res = 'residual: ' + str(stats[0][0][stats[1][1]]) #res
    
    sans.pars.update(parsUpdate[stats[1][1]])
    minParams = sans.pars
    out = [minParams, minPars, stats, res, chi2]
    
    fileType = '.dat'
    fileType2 = '.txt'
    
    save = input('would you like to save the data? Enter 1 for yes [y]: ')
    if save == '1' or save == '':
        location = '../1D_simFits/ReducedChi2_fits/' 
        sans1d.write_3_column(location + sans.description + '_sim' + fileType, sans.expData.x, minSim)
        sans1d.saveStats(location + sans.description + '_stats' + fileType2, minParams=minParams, minPars=minPars, chi2=chi2, res=res)
    
    return out

out = sansFit(options, save)