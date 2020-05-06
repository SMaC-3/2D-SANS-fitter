#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 24 Apr 2020

@author: jkin0004

"""

import numpy as np
from scipy.interpolate import griddata
import csv
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm
from mpl_toolkits import mplot3d
import itertools
from os import path

import sasmodels.compare
import sasmodels.core
import sasmodels.direct_model
import sasmodels.data

import AnnularSectorExtraction_V2 as ansect


# =============================================================================
# Defining grid class
# =============================================================================

class sans2d:  # Necessary? Making another class to be used by two other classes? Changed 2 classes into 3!

    def __init__(self):

        # General initial parameters for performing simulation
        # =============================================================================

        self.qmax = 0.19
        self.qmin = 0.007
        self.nq = 200
        self.skipRows = 4
        self.bandVal = 1

        self.pars = dict(
            scale=1,
            background=0.32,
            sld=0.1,
            sld_solvent=6.4,
            radius=20.0,
            radius_pd=.12,
            radius_pd_n=35,
            #            radius_cap=21,
            length=110.0,
            length_pd=0.16,
            length_pd_n=35,
            theta=90,
            theta_pd=0,
            theta_pd_n=35,
            theta_pd_type='gaussian',
            phi=0,
            phi_pd=18,
            phi_pd_n=35,
            phi_pd_type='gaussian',
            radius_effective_mode=0,
            radius_effective=35,
            volfraction=0.2,
            charge=35,
            temperature=298.0,
            concentration_salt=0.38,
            dielectconst=80.2)

        self.staticPars = dict(
            scale=1.1,
            background=0.4,
            sld=0.1,
            sld_solvent=6.4,
            radius=19.33,
            radius_pd=.12,
            radius_pd_n=35,
            #            radius_cap=21,
            length=88.0,
            length_pd=0,
            length_pd_n=35,
            theta=90,
            theta_pd=90,
            theta_pd_n=35,
            theta_pd_type='uniform',
            phi=0,
            phi_pd=90,
            phi_pd_n=35,
            phi_pd_type='uniform',
            radius_effective_mode=0,
            radius_effective=34,
            volfraction=0.2,
            charge=33,
            temperature=298.0,
            concentration_salt=0.38,
            dielectconst=80.2)

        self.q = np.linspace(-self.qmax, self.qmax, self.nq)
        self.xxq, self.yyq = np.meshgrid(self.q, self.q)

    def getData(self, indexSelect):

        # General method
        # retrieve experimental data using index, mask it, then define sasmodels data object
        # =============================================================================

        csv_filename = '../ExpDataFiles.csv'
        fieldnames = ['index', 'filename', 'sample', 'shear']

        parameters = files_list_reduce(csv_filename, fieldnames)
        fileName, sample, shear = files_to_reduce(parameters, indexSelect)

        self.setName = '../2D_data_sub/' + fileName[0]

        self.expData_raw = np.loadtxt(self.setName, delimiter="  ", skiprows=self.skipRows)
        bsq = np.sqrt(self.expData_raw[:, 0]**2 + self.expData_raw[:, 1]**2)
        self.beamStop = (bsq < self.qmin)
#        self.beamStop = np.logical_or(bsq < self.qmin, bsq>0.15)
        self.expData_bs = self.expData_raw[~self.beamStop]


#        qx_neg_lower = -0.019
#        qx_neg_upper = -0.012
#        qx_pos_lower = 0.021
#        qx_pos_upper = 0.029
#
#        qy_lower = -0.014
#        qy_upper = -0.011
#
#        qx_lower = -0.058
#        qx_upper = -0.051
#
#        qx_neg = np.logical_and(self.expData_bs[:,0] < qx_neg_upper, self.expData_bs[:,0] > qx_neg_lower)
##        qx_neg = (self.expData_bs[:,0]< qx_neg_upper)
#        data_qx_neg = self.expData_bs[~qx_neg]
#
#        qx_pos = np.logical_and(data_qx_neg[:,0] < qx_pos_upper, data_qx_neg[:,0] > qx_pos_lower)
#        data_qx_pos = data_qx_neg[~qx_pos]
#
#        qy_qxRange = np.logical_and(data_qx_pos[:,0] < qx_pos_lower, data_qx_pos[:,0] > qx_neg_upper)
#        qyRange = np.logical_and(data_qx_pos[:,1] < qy_upper, data_qx_pos[:,1] > qy_lower)
#        qy_mask = np.logical_and(qy_qxRange == True, qyRange == True)
#
#        data_qy_mask = data_qx_pos[~qy_mask]
#
#        qx_mask = np.logical_and(data_qy_mask[:,0] < qx_upper, data_qy_mask[:,0] > qx_lower)
#
#        self.expData_masked = data_qy_mask[~qx_mask]

        # Removing detector shadow

        qx_neg_lower = -0.019
        qx_neg_upper = -0.012
        qx_pos_lower = 0.021
        qx_pos_upper = 0.029

        qy_lower = -0.014
        qy_upper = -0.011

        qx_neg = np.logical_and(self.expData_bs[:, 0] <
                                qx_neg_upper, self.expData_bs[:, 0] > qx_neg_lower)
#        qx_neg = (self.expData_bs[:,0]< qx_neg_upper)
        data_qx_neg = self.expData_bs[~qx_neg]

        qx_pos = np.logical_and(data_qx_neg[:, 0] < qx_pos_upper, data_qx_neg[:, 0] > qx_pos_lower)
        data_qx_pos = data_qx_neg[~qx_pos]

        qy_qxRange = np.logical_and(
            data_qx_pos[:, 0] < qx_pos_lower, data_qx_pos[:, 0] > qx_neg_upper)
        qyRange = np.logical_and(data_qx_pos[:, 1] < qy_upper, data_qx_pos[:, 1] > qy_lower)
        qy_mask = np.logical_and(qy_qxRange == True, qyRange == True)

        self.expData_masked = data_qx_pos[~qy_mask]

        # End removing detector shadow

        self.expData_sort = np.array(sorted(self.expData_masked, key=lambda col: (col[1], col[0])))

        self.qx_unique = np.unique(np.array(self.expData_sort[:, 0]))
        self.qy_unique = np.unique(np.array(self.expData_sort[:, 1]))

#        self.expData_sort = abs(self.expData_sort)

        # testing resolution
        self.expData = sasmodels.data.Data2D(x=self.expData_sort[:, 0], dx=0.1*abs(self.expData_sort[:, 0]), y=self.expData_sort[:, 1], dy=0.1*abs(
            self.expData_sort[:, 1]), z=self.expData_sort[:, 2], dz=self.expData_sort[:, 3])
        self.expData.sample = sample
        self.expData.shear = shear

        return

    def interpData(self, data, dataSet=None):

        # General method
        # Interpolate dataSet (2D sasmodels object) onto grid defined using initialised
        # values. Mask out <qmin on grid and I
        # =============================================================================

        #        zzq = griddata(dataSet[:,0:2], dataSet[:,2], (self.xxq,self.yyq), method = 'linear')
        if dataSet != None:
            zzq1 = griddata((dataSet.qx_data, dataSet.qy_data), dataSet.data,
                            (self.xxq, self.yyq), method='linear')
        else:
            zzq1 = griddata((self.expData.qx_data, self.expData.qy_data),
                            data, (self.xxq, self.yyq), method='linear')

        bsq = np.sqrt(self.xxq**2 + self.yyq**2)
        mask = (bsq < self.qmin)

#        zzq = np.full_like(self.xxq, np.nan)
        zzq = np.full_like(self.xxq, 0)
        zzq[~mask] = zzq1[~mask]

        qx_grid = []
        qy_grid = []
        I_grid = []

        for vals in self.xxq:
            for valss in vals:
                qx_grid.append(valss)

        for vals in self.yyq:
            for valss in vals:
                qy_grid.append(valss)

        for vals in zzq:
            for valss in vals:
                I_grid.append(valss)

        qx_grid = np.array(qx_grid)
        qy_grid = np.array(qy_grid)
        I_grid = np.array(I_grid)

        bsq = np.sqrt(qx_grid**2 + qy_grid**2)
        beamStop = (bsq < self.qmin)
        I_grid = I_grid[~beamStop]
        qx_grid = qx_grid[~beamStop]
        qy_grid = qy_grid[~beamStop]

        interp = sasmodels.data.Data2D(x=qx_grid, y=qy_grid, z=I_grid)
        interp.err_data = np.array(list(itertools.repeat(10, len(qx_grid))))

        return interp, zzq
#        return interp

    def getSim(self, indexSelect):

        # General method
        # Sim equivalent of getData. Retrieves simulated data using index, masks it, then defines sasmodels object
        # =============================================================================

        csv_filename = '../SimDataFiles.csv'
        fieldnames = ['index', 'filename', 'sample', 'shear']

        parameters = files_list_reduce(csv_filename, fieldnames)
        fileName, sample, shear = files_to_reduce(parameters, indexSelect)

        self.simName = '../2D_simFits/' + fileName[0]

        self.simImport_raw = np.loadtxt(self.simName, delimiter="  ", skiprows=self.skipRows)
        bsq = np.sqrt(self.simImport_raw[:, 0]**2 + self.simImport_raw[:, 1]**2)
#        self.beamStop = (bsq < self.qmin)
        self.beamStop = np.logical_or(bsq < self.qmin, bsq > 0.15)
        self.simImport_masked = self.simImport_raw[~self.beamStop]

        self.simImport_sort = np.array(
            sorted(self.simImport_masked, key=lambda col: (col[1], col[0])))

        self.simImport = sasmodels.data.Data2D(
            x=self.simImport_sort[:, 0], y=self.simImport_sort[:, 1], z=self.simImport_sort[:, 2])

        return

    def buildSimData(self):

        # General method
        # Experimental data is an imcomplete grid. This method is used to define regions to mask according to the (qx, qy) found
        # from the combination of all unique qx and qy that are not present in the experimental data.
        # A sasmodels data object is defined and this is passed into the makeCalc method
        # =============================================================================

        # Needed for calculator? using q from class definition as want all calcs to be based on the same grid
        self.simData = sasmodels.data.empty_data2D(self.qx_unique, self.qy_unique, resolution=0.1)

        # from here
        # ADD IF STATEMENT TO DETERMINE IF EXP DATA IS PRESENT

        simSet = set(itertools.product(self.qx_unique, self.qy_unique))
        expSet = {(vals[0], vals[1]) for vals in self.expData_sort}

        setDif = simSet - expSet  # (x,y) pairs in simulation not in experiment
        # print(setDif)
        arrayDif = np.array(list(setDif))
        # print(arrayDif)
        ranI = 1000
        # (x,y) pairs in simulation not in experiment with constant z
        self.arrayDif_z = np.insert(arrayDif, 2, ranI, axis=1)
        # print(arrayDif_z)
        exp_data = self.expData_sort[:, 0:3]
        # print(exp_data)
        self.expData_fill = np.vstack((exp_data, self.arrayDif_z))
        self.expData_fill_sort = np.array(
            sorted(self.expData_fill, key=lambda col: (col[1], col[0])))

        # print(exp_dataFill)
        exp_dataFillSas = sasmodels.data.Data2D(
            self.expData_fill_sort[:, 0], self.expData_fill_sort[:, 1], self.expData_fill_sort[:, 2])

        self.simData.mask = np.logical_or(
            exp_dataFillSas.data == ranI, exp_dataFillSas.q_data < self.qmin)

        self.simData.accuracy = 'low'

        # til here could be undefr the getData method

        return

    def buildSasGrid(self):

        # General method
        # Build an empty 2D grid for doing simulations. Can be passed into makeCalc method
        # =============================================================================

        self.q = np.linspace(-self.qmax, self.qmax, self.nq)
        self.sasGrid = sasmodels.data.empty_data2D(self.q)

        return

    def makeCalc(self, dataSet):

        # General method
        #
        # =============================================================================

        #        cyl = sasmodels.core.load_model_info('cylinder')
        #hs = sasmodels.core.load_model_info('hardsphere')
        #        cylhs = sasmodels.core.load_model_info('cylinder@hardsphere')
        cylhmsa = sasmodels.core.load_model_info('cylinder@hayter_msa')
#        capcylhmsa = sasmodels.core.load_model_info('capped_cylinder@hayter_msa')
#        fcylhs = sasmodels.core.load_model_info('flexible_cylinder@hardsphere')

        # model = sasmodels.core.build_model(cylhs, platform = 'dll') #Build using c version instead of python. Avoids pyopencl
        # Build using c version instead of python. Avoids pyopencl
        model = sasmodels.core.build_model(cylhmsa, platform='dll')

        self.calculator = sasmodels.direct_model.DirectModel(dataSet, model)

        return

    def maskq(self, scat):

        # General method
        #
        # =============================================================================

        qx_mask = self.simData.qx_data[~self.simData.mask]
        qy_mask = self.simData.qy_data[~self.simData.mask]

        self.simMask = sasmodels.data.Data2D(
            x=qx_mask, y=qy_mask, z=scat, dz=np.array(list(itertools.repeat(10, len(qx_mask)))))

        return

    def sasPlot(self, data, sim=None, resid=None):

        # General method
        #
        # =============================================================================

        plt.figure()
        sasmodels.data.plot_theory(data, sim, resid, use_data=True, view='log')
        # plt.show()
#        plt.savefig('../2d_simFits/fitPlots/' + self.expData.sample[0][0:-1] + self.expData.shear[0][0:4] +'2D.png')

        return

    def interpPlot(self):

        # General method
        #
        # =============================================================================

        self.expInt, a = self.interpData(self.expData)
        self.simInt, b = self.interpData(self.simMask)
        self.residual = abs(self.expInt.data - self.simInt.data)

        plt.figure()
        self.sasPlot(self.expInt, sim=self.simInt.data, resid=self.residual)

        return

    def sectorCompPlot(self, exp, sim, pltErr=0):

        #        q, binave, errave = ansect.sector(self., centre, width, '0', describer = None)

        expVertq, expVertI, expVerterr = ansect.sector(exp, 0, np.pi/20, '0', describer=None)
        simVertq, simVertI, simVerterr = ansect.sector(sim, 0, np.pi/20, '0', describer=None)
#        print(np.nansum( (( (expVertI-simVertI) / expVerterr)**2 ) ))

        expHorizq, expHorizI, expHorizerr = ansect.sector(
            exp, np.pi/2, np.pi/20, '0', describer=None)
        simHorizq, simHorizI, simHorizerr = ansect.sector(
            sim, np.pi/2, np.pi/20, '0', describer=None)
#        print(np.nansum( ( (expHorizI-simHorizI) / expHorizerr)**2 )/len(expHorizI))

        expq, expI, experr = ansect.sector(exp, np.pi/2, np.pi, '0', describer=None)
        simq, simI, simerr = ansect.sector(sim, np.pi/2, np.pi, '0', describer=None)
#        print(np.nansum( ( (expI-simI) / experr)**2 )/len(expI))

        if pltErr == 1:
            fig = plt.figure(figsize=[8.5, 3], dpi=150)
            ax1 = fig.add_subplot(1, 3, 1)
            ax1.errorbar(expVertq, expVertI, yerr=expVerterr, marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                         label='dataVert')
            ax1.plot(simVertq, simVertI, marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                     label='simVert')
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.legend()
            ax1.minorticks_off()

            ax2 = fig.add_subplot(1, 3, 2)
            ax2.errorbar(expHorizq, expHorizI, yerr=expHorizerr, marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                         label='dataHoriz')
            ax2.plot(simHorizq, simHorizI, marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                     label='simHoriz')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.legend()
            ax2.minorticks_off()

            ax3 = fig.add_subplot(1, 3, 3)
            ax3.plot(expq, expI, marker='o', markersize=2, markerfacecolor=[
                     0, 0, 0], linestyle='', label='dataRad')
            ax3.plot(simq, simI, marker='o', markersize=2, markerfacecolor=[
                     0, 0, 0], linestyle='', label='simRad')
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.legend()
            ax3.minorticks_off()

        else:
            fig = plt.figure(figsize=[8.5, 3], dpi=100)
            ax1 = fig.add_subplot(1, 3, 1)
            ax1.plot(expVertq, expVertI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='dataVert')
            ax1.plot(simVertq, simVertI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='simVert')
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.legend()
            ax1.minorticks_off()
#            plt.xticks(fontsize=8)

            ax2 = fig.add_subplot(1, 3, 2)
            ax2.plot(expHorizq, expHorizI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='dataHoriz')
            ax2.plot(simHorizq, simHorizI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='simHoriz')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.legend()
            ax2.minorticks_off()
#            plt.xticks(fontsize=8)

            ax3 = fig.add_subplot(1, 3, 3)
            ax3.plot(expq, expI, marker='', markersize=2, markerfacecolor=[
                     0, 0, 0], linestyle='-', label='dataRad')
            ax3.plot(simq, simI, marker='', markersize=2, markerfacecolor=[
                     0, 0, 0], linestyle='-', label='simRad')
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.legend()
            ax3.minorticks_off()
            plt.show()
#            plt.xticks(fontsize=8)
#            plt.savefig('../2d_simFits/fitPlots/' + self.expData.sample[0][0:-1] + self.expData.shear[0][0:4] +'sectors.png')

        return

    def annularCompPlot(self, exp, sim, pltErr=0):

        #        q, binave, errave = ansect.sector(self., centre, width, '0', describer = None)

        expVertq, expVertI = ansect.annular(exp, radius=0.07, thx=0.01)
        simVertq, simVertI = ansect.annular(sim, radius=0.07, thx=0.01)

#        if pltErr == 0:
        fig = plt.figure()
        ax = fig.add_axes([1, 1, 1, 1])

        ax.scatter(expVertq, expVertI, marker='o', label='data')

        ax.scatter(simVertq, simVertI, marker='o', label='sim')
#            ax1.set_xscale('log')
#            ax1.set_yscale('log')
        ax.legend()
        ax.minorticks_off()
#        plt.savefig('../2d_simFits/fitPlots/' + self.expData.sample[0][0:-1] + self.expData.shear[0][0:4] +'annulus.png', dpi=150, orientation='landscape')

    def surfacePlot(self, zzq):

        plt.rcParams["figure.figsize"] = 12.8, 9.6
        # Normalize the colors based on Z value
        norm = plt.Normalize(zzq.min(), zzq.max())
        colors = cm.jet(norm(zzq))
        ax = plt.axes(projection='3d')
        surf = ax.plot_surface(self.xxq, self.yyq, zzq, facecolors=colors, shade=False)
        surf.set_facecolor((0, 0, 0, 0))

        return

    def objFunc(self, p_guess, p_list):

        # Least--squares method
        # Defines an objective funciton to minimize. This is a reduced chi squared statistic
        # =============================================================================

        #        print('objFunc called')
        print(dict(zip(p_list, p_guess)))
        dof = len(self.expData.q_data) - len(p_guess)
#        print(dof)

#        bandName = p_list[0]
#        bandVal = p_guess[0]

        omni_pars = ['scale', 'background', 'radius', 'radius_pd',
                     'radius_effective', 'volfraction']

        if 'bandVal' not in p_list and self.bandVal == 1:  # don't need to update static
            #            print('option 1')
            self.pars.update(dict(zip(p_list, p_guess)))
            sim = self.calculator(**self.pars)

        elif any(item in omni_pars for item in p_list):  # need to update static
            #            print('option 2')
            if 'bandVal' in p_list:  # fitting over bandval, need to remove from dicitonary before updating pars

                self.pars.update(dict(zip(p_list, p_guess)))
                del self.pars['bandVal']

                hold = {items: p_guess[idx] for idx, items in enumerate(p_list)}
                staticDict = {}
                for keys, vals in hold.items():
                    if keys in omni_pars:
                        #                        staticDict.update(dict(items = p_guess[idx]))
                        staticDict.update({keys: vals})
                    else:
                        continue

                self.staticPars.update(staticDict)
#                print(self.staticPars)

#                simShear = self.calculator(**self.pars)
#                simStatic = self.calculator(**self.staticPars)

                sim = p_guess[p_list.index('bandVal')]*self.calculator(**self.pars) + \
                    (1 - p_guess[p_list.index('bandVal')])*self.calculator(**self.staticPars)

            elif 'bandVal' not in p_list:
                #                print('option 3')
                self.pars.update(dict(zip(p_list, p_guess)))

                hold = {items: p_guess[idx] for idx, items in enumerate(p_list)}
                staticDict = {}
                for keys, vals in hold.items():
                    if keys in omni_pars:
                        #                        staticDict.update(dict(items = p_guess[idx]))
                        staticDict.update({keys: vals})
                    else:
                        continue

                self.staticPars.update(staticDict)
#                print(self.staticPars)

#                simShear = self.calculator(**self.pars)
#                simStatic = self.calculator(**self.staticPars)

                sim = self.bandVal*self.calculator(**self.pars) + \
                    (1 - self.bandVal)*self.calculator(**self.staticPars)

        else:  # only fitting pars present, no need to update static
            #            print('option 4')
            if 'bandVal' in p_list:

                self.pars.update(dict(zip(p_list, p_guess)))
#                simShear = self.calculator(**self.pars)

                sim = p_guess[p_list.index('bandVal')]*self.calculator(**self.pars) + \
                    (1 - p_guess[p_list.index('bandVal')])*self.simStatic

            elif 'bandVal' not in p_list:

                self.pars.update(dict(zip(p_list, p_guess)))
#                simShear = self.calculator(**self.pars)

                sim = self.bandVal*self.calculator(**self.pars) + (1 - self.bandVal)*self.simStatic

#
#        self.pars.update(dict(zip(p_list[1:], p_guess[1:])))
##        self.pars.update(dict(zip(p_list, p_guess)))
#        simShear = self.calculator(**self.pars)
#
#
#        sim = bandVal*simShear + (1 - bandVal)*self.simStatic

        err = ((((self.expData.data - sim)**2 / self.expData.err_data**2)))/dof
#        print(np.nansum((  ((  (self.expData.data - sim)**2 / self.expData.err_data**2)))))
        print(np.nansum(err))

        return err

# =============================================================================
# End defintion of classes
# =============================================================================


#############################################################################


def buildPars(listInfo):

    #
    # =============================================================================

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


def fitInput(fitChoose):

    #
    # =============================================================================

    p_list = []
    p_guess = []
    p_bounds = ([], [])
    p_num = []

    for keys, values in fitChoose.items():
        if values[0] == 1:
            p_list.append(keys)
            p_guess.append(values[1])
            p_bounds[0].append(values[2])
            p_bounds[1].append(values[3])
            p_num.append(values[4])
        else:
            continue
    return p_list, p_guess, p_bounds, p_num


#############################################################################
def simCalc(parsUpdate, sans):

    #
    # =============================================================================

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

    #
    # =============================================================================

    dof = len(sans.expData.q_data) - len(listInfo)

    res = [sum((abs((sans.expData_sort[:, 2] - sims))/sims) / dof) for sims in scatloop]
    #res = [sum((abs((sans.expData_sort[:,2] - sims))/(sans.expData_sort[:,2]*0.5+sims*0.5)) / dof) for sims in scatloop]
    chi2 = [sum((((sans.expData_sort[:, 2] - sims)**2)/sims) / dof) for sims in scatloop]

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

    minPars = parsUpdate[minResI]

    return listInfo2, stats, minPars

#############################################################################


def getFilenames(csvFilename):

    #
    # =============================================================================

    files = []

    with open(csvFilename) as csv_file:
        reader = csv.DictReader(csv_file, fieldnames=['index', 'filename'])
        iterRows = iter(reader)
        next(iterRows)
        for row in iterRows:
            files.append(row['filename'])
    return files

#############################################################################
# To replace getFilenames?


def files_list_reduce(filename, fieldnames):
    """ Creat array of input reduction settings """
#
# =============================================================================

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

#
# =============================================================================

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
#
# =============================================================================

    files_to_reduce = []
    sample = []
    shear = []

    if len(evaluate_files) == 0:
        files_to_reduce.extend(parameters)
    else:
        # call function for retrieve the IDs list
        evaluate_files_l = evaluate_files_list(evaluate_files)
        for parameter in parameters:
            if int(parameter['index']) in evaluate_files_l:
                files_to_reduce.append(parameter['filename'])
                sample.append(parameter['sample'])
                shear.append(parameter['shear'])

    return files_to_reduce, sample, shear

#############################################################################


def plot(sim, zzq, bounds):

    # DEPRICIATED

    plt.figure(figsize=[6.3, 4.08], dpi=200)
    plt.pcolormesh(sim.xxq, sim.yyq, zzq, cmap='jet',
                   vmin=bounds[0], vmax=bounds[1], norm=matplotlib.colors.LogNorm())
    plt.title('')
    # set the limits of the plot to the limits of the data
    plt.axis([sim.xxq.min(), sim.xxq.max(), sim.yyq.min(), sim.yyq.max()])
    plt.colorbar()
    #plt.pcolormesh(sim.xxq, sim.yyq, zMask[0], vmin = vmin, vmax=vmax, norm=matplotlib.colors.LogNorm())

#############################################################################


def save(sans, describer, minParams, minPars, stats, location, fitInfo):

    #
    # =============================================================================

    #location = '../2D_simFits/ReducedChi2_fits/20wt_CAPB_SLES_2wt_NaCl/'
    describer = describer

    for idx, char in enumerate(sans.expData.shear[0]):
        if char != ' ':
            continue
        else:
            shearIdx = idx
#            print(char)
            break

    shear = sans.expData.shear[0][0:shearIdx]

    name = sans.expData.sample[0] + '_' + shear + 'ps'
    post1 = '_sim'
    type1 = '.dat'

    saveName1 = name + post1 + describer + '_'

    post2 = '_simInfo'
    type2 = '.txt'

    saveName2 = name + post2 + describer + '_'

    output = []

    output.append('qmin = ' + str(sans.qmin))
    output.append('ftol = ' + str(fitInfo[0]))
    output.append('method = ' + str(fitInfo[1]))
    output.append(' ')

    for key, val in minParams.items():
        output.append(str(key) + '=' + str(val) + ',')
    output.append(' ')

    output.append(' static parameters ')
    for key, val in sans.staticPars.items():
        output.append(str(key) + '=' + str(val) + ',')

    output.append(' ')

    output.append('Fitting_performed_over_the_following_parameters:')
    for key in minPars.keys():
        output.append(str(key))

    output.append('Returned_the_following_goodness_of_fit_measures:')
    # output.append(chi2)
    output.append(stats)
#    print(output)

    while path.exists(location) == False:
        print('error: file path does not exist. Please input a valid file path')
        location = input('file path: ')

    if path.exists(location + saveName2 + type2):
        #        versionNum1 = input("Input a version number: ")
        versionNum1 = 'R2'
        with open(location + saveName2 + versionNum1 + type2, 'w') as file:
            for lines in output:
                file.write(lines)
                file.write("\n")
    else:
        #        versionNum1 = input("Input a version number: ")
        versionNum1 = 'R2'
        with open(location + saveName2 + versionNum1 + type2, 'w') as file:
            for lines in output:
                file.write(lines)
                file.write("\n")

    if path.exists(location + saveName1 + type1):
        #        versionNum = input("Input a version number: ")
        write_3_column(location + saveName1 + versionNum1 + type1, sans)

    else:
        write_3_column(location + saveName1 + versionNum1 + type1, sans)

    print('file was saved with filename: ' + saveName1 + versionNum1 + type1)
    return

#############################################################################


def write_3_column(filename, sans):
    #
    # =============================================================================

    with open(filename, 'wt') as fh:
        fh.write("Qx  Qy  I\n")
        for x, y, z in zip(sans.expData.qx_data, sans.expData.qy_data, sans.optimSim):
            fh.write("%g  %g  %g\n" % (x, y, z))

#############################################################################


def sectPlot(sans):

    #
    # =============================================================================

    qxy = np.vstack((sans.simInt.qx_data, sans.simInt.qy_data))
    simData = np.vstack((qxy, sans.simInt.data))
    simData = np.transpose(simData)

    eqxy = np.vstack((sans.expInt.qx_data, sans.expInt.qy_data))
    expData = np.vstack((qxy, sans.expInt.data))
    expData = np.transpose(expData)

    [dataVertBin, dataVertAve] = ansect.sector(expData, 0, np.pi/20)
    [simVertBin, simVertAve] = ansect.sector(simData, 0, np.pi/20)

    [dataHorizBin, dataHorizAve] = ansect.sector(expData, np.pi/2, np.pi/20)
    [simHorizBin, simHorizAve] = ansect.sector(simData, np.pi/2, np.pi/20)

    [dataBin, dataAve] = ansect.sector(expData, np.pi/2, np.pi)
    [simBin, simAve] = ansect.sector(simData, np.pi/2, np.pi)

    plt.figure()
    plt.loglog(dataVertBin, dataVertAve, label='dataVert')
    plt.loglog(simVertBin, simVertAve, label='simVert')
    plt.legend()

    plt.figure()
    plt.loglog(dataHorizBin, dataHorizAve, label='dataHoriz')
    plt.loglog(simHorizBin, simHorizAve, label='simHoriz')
    plt.legend()

    plt.figure()
    plt.loglog(dataBin, dataAve, label='data')
    plt.loglog(simBin, simAve, label='sim')
    plt.legend()

    return
