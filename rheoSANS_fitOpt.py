#!/usr/bin/env python3.
# -*- coding: utf-8 -*-.
"""
Created on Fri 24 Apr 2020

@author: Joshua King

"""
# %% codecell
import numpy as np
from scipy.optimize import least_squares
import rheoSANS_fitOpt_Functions as rsf
import os

# %% codecell
# =============================================================================
# User input
# =============================================================================

# =============================================================================
# Set default parameter dictionaries
# Pars2sim contains parameters to be passed to the cylinder@hayter-msa model for
# the purpose of fitting a high shear, flow aligned band.
# This is the dictionary iterated over in the fitting routine.
# Initial guesses for fitting parameters are taken from this dictionary.
#
# Pars2static contains parameters for a static band. It is recommended that these
# parameters are taken from fitting radially averaged static data with the
# cylinder@hayter-msa model.
#
# "bandVal" is external to these dictionaries are is used to determine the proportions
# of each band.
# =============================================================================

# 0.163

pars2sim = ({'scale':  0.83615,
             'background': 0.43477,
             'sld': -0.4,
             'sld_solvent': 6.3,
             'radius': 19.82,
             'radius_pd': 0.163,
             'radius_pd_n': 10.0,
             'length': 103.44,
             'length_pd': 0.0,
             'length_pd_n': 35.0,
             'theta': 90.0,
             'theta_pd': 90.0,
             'theta_pd_n': 26.0,
             'theta_pd_type': 'uniform',
             'phi': 90.0,
             'phi_pd': 90.0,
             'phi_pd_n': 26.0,
             'phi_pd_type': 'uniform',
             'radius_effective_mode': 1,
             # 'radius_effective': 0,
             'volfraction': 0.216,
             'charge': 30.827,
             'temperature': 298.0,
             'concentration_salt': 0.38,
             'dielectconst': 80.2,
             'radius_pd_type': 'gaussian'})

pars2static = ({'scale':  0.67087,
                'background': 0.35403,
                'sld': -0.4,
                'sld_solvent': 6.3,
                'radius': 19.82,
                'radius_pd': 0.163,
                'radius_pd_n': 10.0,
                'length': 119.76483,
                'length_pd': 0.0,
                'length_pd_n': 35.0,
                'theta': 90.0,
                'theta_pd': 90.0,
                'theta_pd_n': 50.0,
                'theta_pd_type': 'uniform',
                'phi': 0.0,
                'phi_pd': 90.0,
                'phi_pd_n': 50.0,
                'phi_pd_type': 'uniform',
                'radius_effective_mode': 0,
                'radius_effective': 38.194,
                'volfraction': 0.1626,
                'charge': 30.827,
                'temperature': 298.0,
                'concentration_salt': 0.38,
                'dielectconst': 80.2,
                'radius_pd_type': 'gaussian'})

# No radial pd
# pars2sim.update({'radius_pd': 0})
# pars2static.update({'radius_pd': 0})

# Fewer radial pd points
# pars2sim.update({'radius_pd_n': 5.0})
# pars2static.update({'radius_pd_n': 5.0})

# Fewer phi pd points
# pars2sim.update({'phi_pd_n': 20.0})
#
# # Fewer phi & theta pd points in static
# pars2static.update({'phi_pd_n': 20.0})
# pars2static.update({'theta_pd_n': 20.0})

bandVal = 1

# =============================================================================
# Identify experimental data to be used in fitting by referencing index in
# csv file and setting expected concentration and shear rate for error check.
# location will be user specific and will need to modified between users.
# A simple modifcation would be to change this to the current file path.
# =============================================================================

indexSelected = ['71']

conc = '20'  # concentration of sample to be fitted
shear = '0'  # shear rate of sample to be fitted

rsf.input_sample_check(conc, shear, int(indexSelected[0]))
location = rsf.build_save_location(conc, shear)

# =============================================================================
# Select fitting parameters. Initial values taken from pars2sim dictionary.
# =============================================================================

fitChoose = dict(scale=0,
                 background=0,
                 length=0,
                 phi_pd=0,
                 bandVal=0,)

p_list = rsf.fitInput(fitChoose)
p_guess = []
for pars in p_list:
    if pars == 'bandVal':
        p_guess.append(bandVal)
    else:
        p_guess.append(pars2sim[pars])

# =============================================================================
# End user input
# =============================================================================
# =============================================================================
# Defining input
# =============================================================================

options = [indexSelected[0], pars2sim, pars2static, p_list, p_guess,
           location, bandVal]

# =============================================================================
# Defining least squares fitting function
# =============================================================================


def rheoSANS_fitOpt(options, saveOpt):

    # =============================================================================
    # Define items in options
    # =============================================================================

    indexSelect = options[0]
    pars2sim = options[1]  # Update default pars with these
    pars2static = options[2]
    p_list = options[3]
    p_guess = options[4]
    location = options[5]
    bandVal = options[6]

    # =============================================================================
    # Define sans object and perform required methods
    # =============================================================================

    sans = rsf.sans2d()  # Define sans object
    sans.qmin = 0.04  # set qmin
    sans.bandVal = bandVal
    sans.dp = 4  # number of decimal points to which data will be saved
    print('band value = ' + str(sans.bandVal))
    sans.getData(indexSelect)  # get selected experimental data
    # update parameter dictionaries with user input
    sans.pars.update(pars2sim)
    sans.staticPars.update(pars2static)

    # make calculator based on experimental qx, qy. Default resolution = 10%
    sans.makeCalc(sans.expData)

    # =============================================================================
    # Define input parameters for optimiser
    #  =============================================================================

    ftol = 1e-6
    # gtol = 1e-8
    xtol = 1e-4
    max_nfev = 20
    method = 'lm'

    print(sans.expData.sample + sans.expData.shear)
    print(sans.pars)
    print(sans.staticPars)

    # =============================================================================
    # Run optimiser. Define fitting parameters at termination and calculate chi^2_r.
    # Calculate optimal model values and create 2D sasmodels object using this.
    # =============================================================================

    if not p_list:
        optim = sans.objFunc(p_guess=p_guess, p_list=p_list)
        minPars = dict(zip(p_list, p_guess))
        chi2_reduced = ['reduced chi2: ' + str(round(np.sum(optim), sans.dp))]

    else:
        optim = least_squares(sans.objFunc, p_guess, method=method,
                              xtol=xtol, ftol=ftol, args=(p_list,))
        minPars = dict(zip(p_list, optim.x))
        chi2_reduced = ['reduced chi2: ' + str(round(np.sum(optim.fun), sans.dp))]

    sans.pars.update(minPars)
    minPars_complete = sans.pars

    if 'bandVal' in minPars.keys():
        optimSim = minPars['bandVal']*sans.calculator(**sans.pars) + (
            1 - minPars['bandVal'])*sans.calculator(**sans.staticPars)
    else:
        optimSim = sans.bandVal*sans.calculator(**sans.pars) + \
            (1 - sans.bandVal)*sans.calculator(**sans.staticPars)
        minPars_complete.update({'bandVal': sans.bandVal})
    print(minPars)

    sans.optimSim = optimSim
    sans.makeSimObj(optimSim)

    # =============================================================================
    # Calculate statistics of 1D extractions
    # =============================================================================

    chi2_reduced.append(rsf.extract_sector(sans, 0, np.pi/20, 'vertical'))
    chi2_reduced.append(rsf.extract_sector(sans, np.pi/2, np.pi/20, 'horizontal'))
    chi2_reduced.append(rsf.extract_sector(sans, np.pi/2, np.pi, 'radial average'))
    chi2_reduced.append(rsf.extract_annulus(sans, 0.07, 0.01, 'annulus'))

    print(chi2_reduced)

    # =============================================================================
    # Create plots
    # =============================================================================

    # Plot non--interpolated experimental and model data
    sans.sasPlot(data=sans.expData, sim=optimSim)
    # Plot interpolated experimental and model data with residuals
    sans.interpPlot()
    # Plot 1D extractions
    sans.sectorCompPlot(sans.expInt, sans.simInt)
    sans.annularCompPlot(sans.expInt, sans.simInt)

    # =============================================================================
    # Save data
    # =============================================================================

    # os.system('afplay /System/Library/Sounds/Glass.aiff')

    saveOpt = []
    saveOpt = input("would you like to save? enter '1' for yes: ")

    if saveOpt == '1' or saveOpt == '':
        description = input('enter description: ')
        fitInfo = [ftol, method]
        rsf.save(sans, '', minPars_complete, minPars, chi2_reduced, location, fitInfo,
                 description)

    out = [optim, chi2_reduced, sans]

    return out

# =============================================================================
# Run least squares fitting function
# =============================================================================


out = rheoSANS_fitOpt(options, [])

# projPath = '../2D_simFits/ReducedChi2_fits/25wt_CAPB_SLES_2wt_NaCl/25wt_100ps/25wt_100ps_simInfo_preF2.txt'
# pars, pars_static = rsf.loadDict(projPath,[])
# print(pars)
# print(pars)
