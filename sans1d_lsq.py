#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:38:36 2020

@author: jkin0004
"""

import sans1d
import numpy as np

from scipy.optimize import least_squares

# model.parameters.kernel_parameters

# cylinder
pars2sim = (dict(scale=0.8356,
                 background=0.2733,
                 sld=-0.4,
                 sld_solvent=6.3,
                 radius=19.82,
                 radius_pd=0.163,
                 radius_pd_n=35,
                 length=219.5509,
                 length_pd=0,
                 length_pd_n=35,
                 radius_effective_mode=1,
                 radius_effective=35,
                 volfraction=0.0365,
                 charge=30.827,
                 temperature=298.0,
                 concentration_salt=0.38,
                 dielectconst=80.2,
                 ))


indexSelected = '3'
# p_list = ['scale', 'background', 'length']
# p_guess = [0.183, 0.21, 82]

p_list, p_guess = [], []
# p_guess = []

p_bounds = []

options = [pars2sim, p_list, p_guess, p_bounds, indexSelected]


def sansFit(options):

    pars2sim = options[0]  # Update default pars with these
    p_list = options[1]
    p_guess = options[2]
    p_bounds = options[3]
    indexSelect = options[4]

    sans = sans1d.Sans1d()
    sans.qmin = 0.0
    sans.getData(indexSelect)
    sans.pars.update(pars2sim)
    # sans.buildSimData()
    sans.makeCalc(sans.expData)
    sans.dp = 4

    if not p_list:
        optimSim = sans.calculator(**sans.pars)
        out = []
        pre_chi2 = sum(((sans.expData.y - optimSim)**2/sans.expData.dy**2)/len(sans.expData.y))
        chi2 = str(round(pre_chi2, sans.dp))
        minPars = {}
    else:
        optim = least_squares(sans.objFunc, p_guess, ftol=1e-15,
                              xtol=1e-15, method='lm', args=(p_list,))
        print(optim.x)
        out = optim
        minPars = dict(zip(p_list, optim.x))
        sans.pars.update(dict(zip(p_list, optim.x)))
        optimSim = sans.calculator(**sans.pars)
        chi2 = str(round(np.nansum(optim.fun), sans.dp))

    RMSE = (sum((sans.expData.y - optimSim)**2)
            / len(sans.expData.y))**0.5

    sans.sasPlot(data=sans.expData, sim=optimSim)

    fileType = '.dat'
    fileType2 = '.txt'

    save = input('Would you like to save the data? Input 1 for yes [y]: ')

    if save == '1' or save == '':
        version = input('Enter a file descrition: ')
        location = '../1D_simFits/ReducedChi2_fits/'
        sans1d.write_3_column(location + sans.description + '_' + version + '_sim' +
                              fileType, sans.expData.x, optimSim)
        sans1d.saveStats(sans, location + sans.description + '_' + version + '_stats' + fileType2, minParams=sans.pars,
                         minPars=minPars, chi2='chi2_r: ' + chi2, res='RMSE: ' + str(round(RMSE, sans.dp)))

    return out


out = sansFit(options)
