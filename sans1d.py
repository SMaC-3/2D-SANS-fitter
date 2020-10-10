#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 11:53:50 2020

@author: jkin0004
"""

import numpy as np
import sasmodels.compare
import sasmodels.core
import sasmodels.direct_model
import sasmodels.data
import csv
import itertools


class Sans1d:

    def __init__(self):
        self.qmax = 0.19
        self.qmin = 0.006
        self.nq = 50
        self.skipRows = 1

        self.pars = dict(
            scale=1,
            background=0.32,
            sld=0.1,
            sld_solvent=6.4,
            radius=20.0,
            radius_pd=.12,
            radius_pd_n=35,
            #            radius_cap = 20,
            length=110.0,
            length_pd=0,
            length_pd_n=35,
            radius_effective_mode=0,
            radius_effective=35,
            volfraction=0.2,
            charge=35,
            temperature=298.0,
            concentration_salt=0.38,
            dielectconst=80.2)

    def getData(self, indexSelect):

        csv_filename = '../ExpDataFiles_1D_py.csv'
        fieldnames = ['index', 'filename', 'description']

        parameters = files_list_reduce(csv_filename, fieldnames)
        fileName, description = files_to_reduce(parameters, indexSelect)

        self.description = description[0]

        self.setName = '../2D_annular_sector_extraction/py_sect_radAve_exp/' + fileName[0]

        self.expData_raw = np.loadtxt(self.setName, delimiter="  ", skiprows=self.skipRows)

        self.beamStop = (self.expData_raw[:, 0] < self.qmin)
        self.expData_bs = self.expData_raw[~self.beamStop]

        self.expData = sasmodels.data.Data1D(
            x=self.expData_bs[:, 0], y=self.expData_bs[:, 1], dy=self.expData_bs[:, 2])
#        self.expData.err_data = np.array(list(itertools.repeat(10, len(self.expData.y))))
        self.expData.description = description

        return

    def sasPlot(self, data, sim=None):

        sasmodels.data.plot_theory(data, sim, None, use_data=True, view='log')

        return

    def buildSimData(self):

        self.simData = sasmodels.data.empty_data1D(self.expData.x)
        self.simData.accuracy = 'low'
        self.simData.mask = self.simData.x < self.qmin

        return

    def buildSas(self):

        self.q = np.linspace(0, self.qmax, self.nq)
        self.sasGrid = sasmodels.data.empty_data1D(self.q)

        return

    def makeCalc(self, dataSet):

        #cyl = sasmodels.core.load_model_info('cylinder')
        #hs = sasmodels.core.load_model_info('hardsphere')
        #cylhs = sasmodels.core.load_model_info('cylinder@hardsphere')
        capcylhmsa = sasmodels.core.load_model_info('capped_cylinder@hayter_msa')
        cylhmsa = sasmodels.core.load_model_info('cylinder@hayter_msa')

        # model = sasmodels.core.build_model(cylhs, platform = 'dll') #Build using c version instead of python. Avoids pyopencl
        # Build using c version instead of python. Avoids pyopencl
        model = sasmodels.core.build_model(cylhmsa, platform='dll')

        self.calculator = sasmodels.direct_model.DirectModel(dataSet, model)

        return

    def objFunc(self, p_guess, p_list):

        dof = len(self.expData.x) - len(p_list)

        self.pars.update(dict(zip(p_list, p_guess)))
        sim = self.calculator(**self.pars)
        err = ((self.expData.y - sim)**2/self.expData.dy**2)/dof
        print(np.nansum(err))

        return err

#############################################################################


def buildPars(listInfo):

    # Linspace info: parameter, [min,max,numPts]

    # Build a list of dictionaries containing all combos

    numParas = np.shape(listInfo)[0]

    paras = []
    lists = []

    paraInfo = {}

    i = 0

    while i < numParas:
        paras.append(listInfo[i][0][0])
        (numList, step) = np.linspace(listInfo[i][1][0],
                                      listInfo[i][1][1], listInfo[i][1][2], retstep='true')

#        lists.append( numList[1:-1] )
        lists.append(numList)

        paraInfo.update({listInfo[i][0][0] + '_npts': listInfo[i]
                         [1][2], listInfo[i][0][0] + '_steps': step})

        i = i+1

    paraComb = list(itertools.product(*lists))
    paraRep = list(itertools.repeat(paras, len(paraComb)))

    updatePars = []

    print(paraInfo)

    for i in range(len(paraComb)):
        updatePars.append(dict(zip(paraRep[i], paraComb[i])))

    return updatePars, paraInfo

#############################################################################


def simCalc(parsUpdate, sans):
    scatloop = []

    # go through all fitting combos, update dictionary and calculate scattering for each
    # returns column of intensity vals for each calc

    for i, gridVals in enumerate(parsUpdate):
        sans.pars.update(gridVals)
        scatloop.append(sans.calculator(**sans.pars))
        print(str(i+1) + ' of ' + str(len(parsUpdate)))

    return scatloop

#############################################################################


def errorCalc(listInfo, scatloop, parsUpdate, paraInfo, sans):

    dof = len(sans.expData.x) - len(listInfo)

    res = [sum((abs((sans.expData_bs[:, 1] - sims))/sims) / dof) for sims in scatloop]
    #res = [sum((abs((sans.expData_sort[:,2] - sims))/(sans.expData_sort[:,2]*0.5+sims*0.5)) / dof) for sims in scatloop]
#    chi2 = [sum(  (((sans.expData_bs[:,1] - sims)**2)/sims) / dof ) for sims in scatloop]

    chi2 = [np.nansum(((sans.expData.y - sims)**2 / sans.expData.dy**2))/dof for sims in scatloop]

    #minRes = [np.amin(zzqResSum) for set in range(np.shape(zzqRes)[0])]
    minResI = np.argmin(res)
    minChi2I = np.argmin(chi2)

    stats = [[res, minResI], [chi2, minChi2I]]

    # print(absMinResI)

    # print(parsUpdate[absMinResI])

    listInfo2 = []

    for para, vals in parsUpdate[minResI].items():
        listInfo2.append([[para], [vals - paraInfo[para + '_steps'], vals +
                                   paraInfo[para + '_steps'], paraInfo[para + '_npts']]])

#    minPars = parsUpdate[minResI]
    minPars = parsUpdate[minChi2I]

    return listInfo2, stats, minPars

##############################################################################


def files_list_reduce(filename, fieldnames):
    """ Creat array of input reduction settings """
    parameters = []
    with open(filename) as csv_file:

        reader = csv.DictReader(csv_file, fieldnames=fieldnames)
        iterRows = iter(reader)
        next(iterRows)
        for row in iterRows:
            if row['index'] == '':
                continue
            if row['index'] == 'END':
                break
            parameters.append(row)
    return parameters

##############################################################################


def evaluate_files_list(numbers):
    """ Needed for FilesToReduce, see below """

    expanded = []
    for number in numbers.split(","):
        if "-" in number:
            start, end = number.split("-")
            nrs = range(int(start), int(end) + 1)
            expanded.extend(nrs)
        else:
            expanded.append(int(number))
    return expanded

##############################################################################


def files_to_reduce(parameters, evaluate_files):
    """ Create list of the files to reduce """

    files_to_reduce = []
    description = []

    if len(evaluate_files) == 0:
        files_to_reduce.extend(parameters)
    else:
        # call function for retrieve the IDs list
        evaluate_files_l = evaluate_files_list(evaluate_files)
        for parameter in parameters:
            if int(parameter['index']) in evaluate_files_l:
                files_to_reduce.append(parameter['filename'])
                description.append(parameter['description'])

    return files_to_reduce, description

##############################################################################


def write_3_column(filename, q, I):
    with open(filename, 'wt') as fh:
        fh.write("Q I(Q)\n")
        for x, y in zip(q, I):
            fh.write("%g  %g\n" % (x, y))


def saveStats(sans, filename, minParams, minPars, chi2, res):
    output = []

    for key, val in minParams.items():
        if type(val) == str:
            output.append(str(key) + '=' + str(val) + ',')
        else:
            output.append(str(key) + '=' + str(round(val, sans.dp)) + ',')

    # for key, val in minParams.items():
    #     if
    #     output.append(str(key) + '=' + str(val) +',')
#    for key, val in minParams.items():
#        output.append(str(key) + ' ' + str(val))

    output.append('Fitting_performed_over_the_following_parameters:')
    for key in minPars.keys():
        output.append(str(key))

    output.append('Returned_the_following_goodness_of_fit_measures:')
    output.append(chi2)
    output.append(res)

    with open(filename, 'wt') as file:
        for lines in output:
            file.write(lines)
            file.write("\n")
