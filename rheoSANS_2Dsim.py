#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 11:56:42 2020

@author: jkin0004
"""
import rheoSANSFunctions_fitOpt_lsq as rsf
import os
from scipy.optimize import least_squares
import numpy as np

saveOpt = '0'  # set to '1' to save processed data files, set to any 'number string' other than '1' to not save
describer = '20pc_static_1Dpars'

# ['1' (version #), indexSelect, pars2sim, [ [fitParam, [fitMin, fitMax, numPt1]] ], loopNum]

pars2sim = (dict(
            scale=1.25,
            background=0.38,
            sld=.9,
            sld_solvent=6.4,
            radius=19.4,
            radius_pd=0,
            radius_pd_n=35,
            #            radius_cap=21,
            length=125,
            length_pd=0.18,
            length_pd_n=35,
            theta=90,
            theta_pd=20,
            theta_pd_n=35,
            theta_pd_type='gaussian',
            phi=0,
            phi_pd=23,
            phi_pd_n=35,
            phi_pd_type='gaussian',
            radius_effective_mode=0,
            radius_effective=36,
            volfraction=0.215,
            charge=40,
            temperature=298.0,
            concentration_salt=0.38,
            dielectconst=80.2))


indexSelected = ['67']


options = [describer, None, pars2sim]


def rheoSANS_fitOpt(options, saveOpt):

    # Define items in options

    describer = options[0]
    indexSelect = options[1]
#    pars2sim = options[2] #Update default pars with these
#    pars2static = options[3]

    # Define sans object and perform required methods

    sans = rsf.sans2d()
    sans.qmin = 0.007
    sans.getData(indexSelect)
    sans.pars.update(pars2sim)
#    sans.staticPars.update(pars2static)
    sans.buildSimData()
    sans.makeCalc(sans.simData)
#    sans.simStatic = sans.calculator(**sans.staticPars)

    optimSim = sans.calculator(**sans.pars)
    sans.maskq(optimSim)  # applies mask to qx and qy (already applied to I(q))

#    Plot experimental and min simulated data, sector analysis
    sans.sasPlot(data=sans.expData, sim=optimSim)
    sans.interpPlot()
    rsf.sectPlot(sans)
#    print(sans.pars)
#    print('band val:' + str(bandVal))

#    os.system('afplay /System/Library/Sounds/Glass.aiff')
#    out = [optimSim, sans]
    dof = len(sans.expData.q_data)
    res = sum((abs((sans.expData_sort[:, 2] - optimSim))/optimSim) / dof)
    chi2 = sum((((sans.expData_sort[:, 2] - optimSim)**2)/optimSim) / dof)
    filename = describer + '.dat'

    rsf.write_3_column(filename, sans)

    return out


out = []
for idx in indexSelected:
    options[1] = idx
    out.append(rheoSANS_fitOpt(options, saveOpt))
