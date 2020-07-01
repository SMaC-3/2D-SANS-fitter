#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 24 Apr 2020

@author: jkin0004

"""
import numpy as np
from scipy.optimize import least_squares
import rheoSANSFunctions_fitOpt_omni as rsf
import AnnularSectorExtraction_V2 as ansect
import os
from lmfit import minimize, Parameters
# from scipy.optimize import basinhopping

# =============================================================================
# User input
# =============================================================================
# CHECKLIST

# SAVE OPT, DESCRIBER AND LOCATION SET?
# INDEX SELECTED CORRECT?
# PARAMETER DICTIONARIES AND FIT LIST CORRECT?
# FTOL APPROPRIATE?

# =============================================================================

# =============================================================================
# Set default parameter dictionaries
# =============================================================================

pars2sim = ({'scale': 0.841778,
             'background': 0.38474,
             'sld': -0.4,
             'sld_solvent': 6.3,
             'radius': 19.82,
             'radius_pd': 0,  # 0.163,
             'radius_pd_n': 35.0,
             'length': 187.37694,
             'length_pd': 0.0,
             'length_pd_n': 35.0,
             'theta': 90.0,
             'theta_pd': 0.0,
             'theta_pd_n': 35.0,
             'theta_pd_type': 'gaussian',
             'phi': 0.0,
             'phi_pd': 17.81125,
             'phi_pd_n': 35.0,
             'phi_pd_type': 'gaussian',
             'radius_effective_mode': 0.0,
             'radius_effective': 33.775,
             'volfraction': 0.268,
             'charge': 30.827,
             'temperature': 298.0,
             'concentration_salt': 0.38,
             'dielectconst': 80.2,
             'radius_pd_type': 'gaussian', })

pars2static = ({'scale': 0.841778,
                'background': 0.38474,
                'sld': -0.4,
                'sld_solvent': 6.3,
                'radius': 19.82,
                'radius_pd': 0,  # 0.163,
                'radius_pd_n': 35.0,
                'length': 91.472,
                'length_pd': 0.0,
                'length_pd_n': 35.0,
                'theta': 90.0,
                'theta_pd': 90.0,
                'theta_pd_n': 35.0,
                'theta_pd_type': 'uniform',
                'phi': 0.0,
                'phi_pd': 90.0,
                'phi_pd_n': 35.0,
                'phi_pd_type': 'uniform',
                'radius_effective_mode': 0.0,
                'radius_effective': 33.775,
                'volfraction': 0.268,
                'charge': 30.827,
                'temperature': 298.0,
                'concentration_salt': 0.38,
                'dielectconst': 80.2,
                'radius_pd_type': 'gaussian'})

# parameters to leave unchanged in above dicts
# popPars = ['phi_pd_type', 'theta_pd_type', 'phi_pd', 'theta_pd']
popPars = ['phi_pd_type', 'theta_pd_type', 'theta_pd']

# =============================================================================
# Set sample specific variables
# =============================================================================

# 84-95
indexSelected = ['77']

bandVal_man = 0.287175

saveOpt = '0'  # set to '1' to save processed data files, set to any 'number string' other than '1' to not save
describer = ''

conc = '25'
shear_u = '100'  # shear rate of sample being used to update
shear_f = shear_u  # shear rate of sample being fitted
version = 'R1'
use2update = 'n'

rsf.input_sample_check(conc, shear_f, int(indexSelected[0]))
location = rsf.build_save_location(conc, shear_f)

if use2update == 'y':

    sim, static, bandVal_man = rsf.update_dict_from_file(
        conc, shear_u, version, popPars=popPars, bandVal_man=bandVal_man)
    pars2sim.update(sim)
    pars2static.update(static)
    # print(pars2sim)
    print('Dictionaries successfully updated from file, excluding ' +
          ', '.join(popPars))
else:
    print('User defined dictionaries used')

# =============================================================================
# Select fitting parameters & initial values
# =============================================================================

fitChoose = dict(scale=[0,     0.76,          0.9,        1.5,        5],
                 background=[0,     0.328,        0.01,       0.6,        5],
                 radius=[0,     21,         18,         22,         5],
                 radius_pd=[0,     0.1,        0,          0.5,        5],
                 length=[1,     133,        80,         250,        5],
                 length_pd=[0,     0.1,        0,          0.9,        5],
                 phi=[0,     0,          0,          90,         5],
                 phi_pd=[0,     45,         0,          90,         5],
                 theta=[0,     90,         0,          90,         5],
                 theta_pd=[0,     20,          0,          90,         5],
                 radius_effective=[0,     36,         20,         60,         5],
                 volfraction=[0,     0.,     0.02,       0.35,       5],
                 charge=[0,     27,         10,         40,         5],
                 bandVal=[0,     0.57276,        0,          1,          5],)

p_list, p_guess, p_bounds, p_num = rsf.fitInput(fitChoose)
p_bounds = []
p_guess = []
for pars in p_list:
    if pars == 'bandVal':
        p_guess.append(fitChoose[pars][1])
    else:
        p_guess.append(pars2sim[pars])

# =============================================================================
# Defining input
# =============================================================================

options = [describer, None, pars2sim, p_list, p_guess, p_bounds, pars2static,
           location, bandVal_man]

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
    bandVal_man = options[8]

# =============================================================================
    # Define sans object and perform required methods

    sans = rsf.sans2d()  # Define sans object
    sans.qmin = 0.04  # set qmin
    sans.bandVal = bandVal_man
    print('band value = ' + str(sans.bandVal))
    sans.getData(indexSelect)  # get selected experimental data
    sans.pars.update(pars2sim)  # update parameter dictionary with user input
    sans.staticPars.update(pars2static)  # update static parameter dictionary with user input
#    sans.buildSimData()         #build data object for simulation masking low q and incomplete grid

    sans.makeCalc(sans.expData)  # make calculator using masked simulation object

    sans.simStatic = sans.calculator(**sans.staticPars)  # calculate static simulation

    ftol = 1e-6
    # gtol = 1e-8
    xtol = 1e-4
    max_nfev = 20
    method = 'lm'

    sans.dp = 4

    print(sans.expData.sample + sans.expData.shear)
    print(sans.pars)
    print(sans.staticPars)
    if not p_list:
        optim = sans.objFunc(p_guess=p_guess, p_list=p_list)
        # print(np.nansum(optim))
    else:
        optim = least_squares(sans.objFunc, p_guess, method=method,
                              xtol=xtol, ftol=ftol, args=(p_list,))
        # params = Parameters()
        # params.add('length', value=190)
        # optim = minimize(sans.objFunc, params=pars2sim)
#    optim = basinhopping(sans.objFunc, p_guess,minimizer_kwargs={'args': (p_list,)})
# args=(p_list,))

# =============================================================================
    # Calculate statistics

#    print(optim.x)
    out = optim

    if not p_list:
        minPars = dict(zip(p_list, p_guess))
    else:
        minPars = dict(zip(p_list, optim.x))
    sans.pars.update(minPars)
    minParams = sans.pars
    # print(minParams)
#    optimSim = sans.calculator(**sans.pars)
    if 'bandVal' in minPars.keys():
        optimSim = minPars['bandVal']*sans.calculator(**sans.pars) + (
            1 - minPars['bandVal'])*sans.calculator(**sans.staticPars)
    else:
        optimSim = sans.bandVal*sans.calculator(**sans.pars) + \
            (1 - sans.bandVal)*sans.calculator(**sans.staticPars)
        minParams.update({'bandVal': sans.bandVal})
    print(minPars)
    sans.optimSim = optimSim

    if not p_list:
        chi2_reduced = ['reduced chi2: ' + str(round(np.sum(optim), sans.dp))]
    else:
        chi2_reduced = ['reduced chi2: ' + str(round(np.sum(optim.fun), sans.dp))]  # res

    sans.makeSimObj(optimSim)
    sectorSettings = [[0, np.pi/20, 'reduced chi2 of vertical sector: '],
                      [np.pi/2, np.pi/20, 'reduced chi2 of horizontal sector: '],
                      [np.pi/2, np.pi, 'reduced chi2 of radial average: ']]
    # chi2_sects = []

    for ssets in sectorSettings:
        q_exp, I_exp, err_exp = ansect.sector(sans.expData, ssets[0], ssets[1])
        q_sim, I_sim, err_sim = ansect.sector(sans.simData, ssets[0], ssets[1])

        np.nansum((((I_exp-I_sim)/err_exp)**2)/len(q_exp))

        statString = ssets[2] + \
            str(round(np.nansum((((I_exp-I_sim)/err_exp)**2)/len(q_exp)), sans.dp))
        chi2_reduced.append(statString)
        # print(len(q_exp))

    q_an_exp, I_an_exp, err_an_exp = ansect.annular(dataSet=sans.expData, radius=0.07, thx=0.01)
    q_an_sim, I_an_sim, err_an_sim = ansect.annular(dataSet=sans.simData, radius=0.07, thx=0.01)

    chi2_reduced.append('reduced chi2 of annulus: ' +
                        str(round(np.nansum((((I_an_exp-I_an_sim)/err_an_exp)**2)/len(q_an_exp)), sans.dp)))
    print(chi2_reduced)
# =============================================================================

    # Plot experimental and min simulated data, sector analysis
    sans.sasPlot(data=sans.expData, sim=optimSim)

    sans.expInt, expZzq = sans.interpData(None, dataSet=sans.expData)
    sans.simInt, simZzq = sans.interpData(optimSim)  # Fix this
    sans.residual = abs(sans.expInt.data - sans.simInt.data)
    resZzq = abs(expZzq - simZzq)

    sans.sasPlot(sans.expInt, sim=sans.simInt.data, resid=sans.residual)

    sans.sectorCompPlot(sans.expInt, sans.simInt)

    sans.annularCompPlot(sans.expInt, sans.simInt)


#    sans.surfacePlot(expZzq)
#    sans.surfacePlot(resZzq)

    os.system('afplay /System/Library/Sounds/Glass.aiff')

    saveOpt = []
    saveOpt = input("would you like to save? enter '1' for yes: ")

    if saveOpt == '1' or saveOpt == '':
        fitInfo = [ftol, method]
        rsf.save(sans, describer, minParams, minPars, chi2_reduced, location, fitInfo)
        # rsf.save(sans,'',sans.pars,dict(zip(p_list, out[0][0].x)),out[0][1],location, [1e-6, 'lm'])

    out = [optim, chi2_reduced, sans]

    return out


# =============================================================================
# Running least squares fitting function
# =============================================================================

out = []
for idx in indexSelected:
    options[1] = idx
    out.append(rheoSANS_fitOpt(options, saveOpt))
