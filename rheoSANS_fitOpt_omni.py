#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 24 Apr 2020

@author: jkin0004

"""
import rheoSANSFunctions_fitOpt_omni as rsf
import os
from scipy.optimize import least_squares
# from scipy.optimize import basinhopping
import numpy as np

# =============================================================================
# User input
# =============================================================================
# CHECKLIST

# SAVE OPT, DESCRIBER AND LOCATION SET?
# INDEX SELECTED CORRECT?
# PARAMETER DICTIONARIES AND FIT LIST CORRECT?
# FTOL APPROPRIATE?

# =============================================================================
print('test')

saveOpt = '0'  # set to '1' to save processed data files, set to any 'number string' other than '1' to not save
describer = ''
location = '../2D_simFits/ReducedChi2_fits/25wt_CAPB_SLES_2wt_NaCl/'


# ['1' (version #), indexSelect, pars2sim, [ [fitParam, [fitMin, fitMax, numPt1]] ], loopNum]

pars2sim = (dict(scale=0.8005,
                 background=0.5142867229475522,
                 sld=-0.4,
                 sld_solvent=6.3,
                 radius=21.350361527275133,
                 radius_pd=0,
                 radius_pd_n=35,
                 length=82.2,
                 length_pd=0,
                 length_pd_n=35,
                 theta=90,
                 theta_pd=0,
                 theta_pd_n=35,
                 theta_pd_type='gaussian',
                 phi=0,
                 phi_pd=20,  # try with larger pd like 360
                 phi_pd_n=35,
                 phi_pd_type='gaussian',
                 radius_effective_mode=0,
                 radius_effective=33.474445715581474,
                 volfraction=0.2683,
                 charge=26.978147871213366,
                 temperature=298.0,
                 concentration_salt=0.38,
                 dielectconst=80.2,))

pars2static = (dict(scale=0.8090981755706477,
                    background=0.5142867229475522,
                    sld=-0.4,
                    sld_solvent=6.3,
                    radius=21.350361527275133,
                    radius_pd=0,
                    radius_pd_n=35,
                    length=82.2,
                    length_pd=0,
                    length_pd_n=35,
                    theta=90,
                    theta_pd=90,
                    theta_pd_n=35,
                    theta_pd_type='uniform',
                    phi=0,
                    phi_pd=90,  # try with larger pd like 360
                    phi_pd_n=35,
                    phi_pd_type='uniform',
                    radius_effective_mode=0,
                    radius_effective=33.474445715581474,
                    volfraction=0.2683,
                    charge=26.978147871213366,
                    temperature=298.0,
                    concentration_salt=0.38,
                    dielectconst=80.2,))

indexSelected = ['61']

# par = [0/1, guess, lower bound, upper bound, num pts]

fitChoose = dict(scale=[0,     0.8005,          0.9,        1.5,        5],
                 background=[0,     0.51,        0.01,       0.6,        5],
                 radius=[0,     20.3,         18,         22,         5],
                 radius_pd=[0,     0.1,        0,          0.5,        5],
                 length=[1,     145,        80,         250,        5],
                 length_pd=[0,     0.1,        0,          0.9,        5],
                 phi=[0,     0,          0,          90,         5],
                 phi_pd=[0,     17,         0,          90,         5],
                 theta=[0,     90,         0,          90,         5],
                 theta_pd=[0,     40,          0,          90,         5],
                 radius_effective=[0,     33.5,         20,         60,         5],
                 volfraction=[0,     0.2157,     0.02,       0.35,       5],
                 charge=[0,     27,         10,         40,         5],
                 bandVal=[0,     0.5,        0,          1,          5],)


#p_list = ['scale']
#p_guess = [0.8]
# p_bounds=[]

p_list, p_guess, p_bounds, p_num = rsf.fitInput(fitChoose)
p_bounds = []

options = [describer, None, pars2sim, p_list, p_guess, p_bounds, pars2static, location]

# =============================================================================
# Defining least squares fitting function
# =============================================================================


def rheoSANS_fitOpt(options, saveOpt):

    # =============================================================================
    # Define items in options

    describer = options[0]
    indexSelect = options[1]
    pars2sim = options[2]  # Update default pars with these
    p_list = options[3]
    p_guess = options[4]
    p_bounds = options[5]
    pars2static = options[6]
    location = options[7]

# =============================================================================
    # Define sans object and perform required methods

    sans = rsf.sans2d()  # Define sans object
    sans.qmin = 0.04  # set qmin
    sans.bandVal = 1
    sans.getData(indexSelect)  # get selected experimental data
    sans.pars.update(pars2sim)  # update parameter dictionary with user input
    sans.staticPars.update(pars2static)  # update static parameter dictionary with user input
#    sans.buildSimData()         #build data object for simulation masking low q and incomplete grid

    sans.makeCalc(sans.expData)  # make calculator using masked simulation object

    sans.simStatic = sans.calculator(**sans.staticPars)  # calculate static simulation

    ftol = 1e6
    method = 'lm'

    optim = least_squares(sans.objFunc, p_guess, method=method, ftol=ftol,
                          args=(p_list,))
#    optim = basinhopping(sans.objFunc, p_guess,minimizer_kwargs={'args': (p_list,)})
# args=(p_list,))

# =============================================================================
    # Calculate statistics

#    print(optim.x)
    out = optim
    minPars = dict(zip(p_list, optim.x))
    sans.pars.update(minPars)
    minParams = sans.pars
    print(minParams)
#    optimSim = sans.calculator(**sans.pars)
    if 'bandVal' in minPars.keys():
        optimSim = minPars['bandVal']*sans.calculator(**sans.pars) + (
            1 - minPars['bandVal'])*sans.calculator(**sans.staticPars)
    else:
        optimSim = sans.bandVal*sans.calculator(**sans.pars) + \
            (1 - sans.bandVal)*sans.calculator(**sans.staticPars)
    sans.optimSim = optimSim

    chi2_reduced = np.sum(optim.fun)

    chi2_reduced = 'reduced chi2: ' + str(np.sum(optim.fun))  # res

# =============================================================================

    # Plot experimental and min simulated data, sector analysis
    sans.sasPlot(data=sans.expData, sim=optimSim)

    sans.expInt, expZzq = sans.interpData(None, dataSet=sans.expData)
    sans.simInt, simZzq = sans.interpData(optimSim)  # Fix this
    sans.residual = abs(sans.expInt.data - sans.simInt.data)
    resZzq = abs(expZzq - simZzq)

    sans.sasPlot(sans.expInt, sim=sans.simInt.data, resid=sans.residual)

    sans.sectorCompPlot(sans.expInt, sans.simInt)
    #
    sans.annularCompPlot(sans.expInt, sans.simInt)


#    sans.surfacePlot(expZzq)
#    sans.surfacePlot(resZzq)

    os.system('afplay /System/Library/Sounds/Glass.aiff')

#    saveOpt = []
#    saveOpt = input("would you like to save? enter '1' for yes: ")

    if saveOpt == '1' or saveOpt == '':
        fitInfo = [ftol, method]
        rsf.save(sans, describer, minParams, minPars, chi2_reduced, location, fitInfo)

    out = [optim, chi2_reduced, sans]

    return out


# =============================================================================
# Running least squares fitting function
# =============================================================================

out = []
for idx in indexSelected:
    options[1] = idx
    out.append(rheoSANS_fitOpt(options, saveOpt))
