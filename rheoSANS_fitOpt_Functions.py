#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 24 Apr 2020

@author: Joshua King

"""

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import csv
import matplotlib.pyplot as plt
# import matplotlib.colors
from matplotlib import cm
import itertools
from os import path
import os

import sasmodels.compare
import sasmodels.core
import sasmodels.direct_model
import sasmodels.data

import annular_sector_extraction as ansect
import datetime


# =============================================================================
# Defining grid class
# =============================================================================

class sans2d:  # Necessary? Making another class to be used by two other classes? Changed 2 classes into 3!
    """Class object for fitting and plotting 2D SANS data."""

    def __init__(self):

        self.qmax = 0.19
        self.qmin = 0.04
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
        """Load experimental scattering data from a given index.

        The index value should refer to a row defined by an 'index' column in a
        csv file with other column headings 'filename', 'sample', and 'shear'.

        Data is loaded from 'filename' corresponding to the given index. The scattering data
        should be formatted as follows: qx, qy, I(q), I_err(q)

        Args:
            indexSelect

        Returns:
            2D sasmodels object
        """
        csv_filename = '../ExpDataFiles.csv'
        fieldnames = ['index', 'filename', 'sample', 'shear']

        parameters = files_list_reduce(csv_filename, fieldnames)
        fileName, sample, shear = files_to_reduce(parameters, indexSelect)

        self.setName = '../2D_data_sub/' + fileName[0]

        self.expData_raw = np.loadtxt(self.setName, delimiter="  ", skiprows=self.skipRows)
        bsq = np.sqrt(self.expData_raw[:, 0]**2 + self.expData_raw[:, 1]**2)
        self.beamStop = (bsq < self.qmin)
        self.expData_bs = self.expData_raw[~self.beamStop]

        # Removing detector shadow

        qx_neg_lower = -0.019
        qx_neg_upper = -0.012
        qx_pos_lower = 0.021
        qx_pos_upper = 0.029

        qy_lower = -0.014
        qy_upper = -0.011

        qx_neg = np.logical_and(self.expData_bs[:, 0] <
                                qx_neg_upper, self.expData_bs[:, 0] > qx_neg_lower)
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
        """By default, this method takes I(q) data and interpolates this using qx, qy
        from previously loaded expData. If dataSet is given instead, then qx, qy data
        from that is used. Masks out q < qmin.

        Args:
            data, dataSet (optional)

        Returns:
            2D sasmodels object, zzq grid of interpolated values.
        """

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
        interp.err_data = np.array(list(itertools.repeat(0, len(qx_grid))))

        return interp, zzq

    def interpData_noMask(self, data, dataSet=None):
        """By default, this method takes I(q) data and interpolates this using qx, qy
        from previously loaded expData. If dataSet is given instead, then qx, qy data
        from that is used. Masks out q < qmin.

        Args:
            data, dataSet (optional)

        Returns:
            2D sasmodels object, zzq grid of interpolated values.
        """

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

        # bsq = np.sqrt(qx_grid**2 + qy_grid**2)
        # beamStop = (bsq < self.qmin)
        # I_grid = I_grid[~beamStop]
        # qx_grid = qx_grid[~beamStop]
        # qy_grid = qy_grid[~beamStop]

        interp = sasmodels.data.Data2D(x=qx_grid, y=qy_grid, z=I_grid)
        interp.err_data = np.array(list(itertools.repeat(0, len(qx_grid))))

        return interp, zzq

    def getSim(self, indexSelect):
        """Load saved simulated scattering data from a given index.

        The index value should refer to a row defined by an 'index' column in a
        csv file with other column headings 'filename', 'sample', and 'shear'.

        Data is loaded from 'filename' corresponding to the given index. The scattering
        data should be formatted as follows: qx, qy, I(q)

        Args:
            indexSelect

        Returns:
            2D sasmodels object
        """
        csv_filename = '../SimDataFiles.csv'
        fieldnames = ['index', 'filename', 'sample', 'shear']

        parameters = files_list_reduce(csv_filename, fieldnames)
        fileName, sample, shear = files_to_reduce(parameters, indexSelect)

        self.simName = '../2D_simFits/ReducedChi2_fits/allConc_bestFits/' + fileName[0]

        self.simImport_raw = np.loadtxt(self.simName, delimiter="  ", skiprows=self.skipRows)
        bsq = np.sqrt(self.simImport_raw[:, 0]**2 + self.simImport_raw[:, 1]**2)
#        self.beamStop = (bsq < self.qmin)
        self.beamStop = bsq < self.qmin
        self.simImport_masked = self.simImport_raw[~self.beamStop]

        self.simImport_sort = np.array(
            sorted(self.simImport_masked, key=lambda col: (col[1], col[0])))

        dz = np.array(list(itertools.repeat(0, len(self.simImport_sort[:, 0]))))
        self.simImport = sasmodels.data.Data2D(
            x=self.simImport_sort[:, 0], y=self.simImport_sort[:, 1],
            z=self.simImport_sort[:, 2], dz=dz)

        return

    def buildSasGrid(self):
        """
        Build an empty 2D grid for doing simulations.

        Can be passed into makeCalc method

        Args:
            None (uses values defined in __init__)

        Returns:
            Empty 2D sasmodels object
        """

        self.q = np.linspace(-self.qmax, self.qmax, self.nq)
        self.sasGrid = sasmodels.data.empty_data2D(self.q)

        return

    def makeCalc(self, dataSet):
        """Build calculator object to perform simulations based on cylinder@hmsa.

        Args:
            dataSet: 2D sasmodels data object

        Returns:
            calculator object
        """

        #cyl = sasmodels.core.load_model_info('cylinder')
        #hs = sasmodels.core.load_model_info('hardsphere')
        #cylhs = sasmodels.core.load_model_info('cylinder@hardsphere')
        cylhmsa = sasmodels.core.load_model_info('cylinder@hayter_msa')

        # Build using c version instead of python. Avoids pyopencl
        model = sasmodels.core.build_model(cylhmsa, platform='dll')
        self.calculator = sasmodels.direct_model.DirectModel(dataSet, model)

        return

    def makeSimObj(self, scat):
        """Make 2D sasmodels object with model calculated I(q).

        This I(q) should have been calculated based on imported experimental qx, qy

        Args:
            scat: model calculated I(q)

        Returns:
            simData 2D sasmodels object
        """

        dz = np.array(list(itertools.repeat(0, len(self.expData_sort[:, 0]))))
        self.simData = sasmodels.data.Data2D(x=self.expData_sort[:, 0],
                                             y=self.expData_sort[:, 1],
                                             z=scat,
                                             dz=dz)
        return

    def sasPlot(self, data, sim=None, resid=None):
        """Plot 2D data, model and/or residuals

        Args:
            data: 2D sasmodels object of experimental data
            sim: 2D grid of I(q) from sasmdoels calculator
            residuals: 2D grid of residuals

        Returns:
            Sasmodels plot/sub plots
        """
        plt.figure()
        sasmodels.data.plot_theory(data, sim, resid, use_data=True, view='log')
        return

    def interpPlot(self):
        """
        Plots interpolated experimental and modelled data with residuals.

        Requires sans2d object to be defined with expData, optimSim attributes.

        Args:
            None

        Returns:
            Sasmodels plot/sub plots using interpolated data
        """
        self.expInt, a = self.interpData(None, self.expData)
        self.simInt, b = self.interpData(self.optimSim)
        self.residual = abs(self.expInt.data - self.simInt.data)

        plt.figure()
        self.sasPlot(self.expInt, sim=self.simInt.data, resid=self.residual)

        return

    def sectorCompPlot(self, exp, sim, pltErr=0):
        """Plots vertical, horizontal and radially averaged sectors from experimental
        and modelled data.

        Args:
            exp: 2D sasmodels object of experimental data
            sim: 2D sasmodels object of model data
            pltErr (optional): default = 0. Set to 1 to plot with errorbars.
            Errorbars are not valid for interpolated data.

        Returns:
            Subplots of vertical, horizontal and radially averaged sectors
        """
        expVertq, expVertI, expVerterr = ansect.sector(
            exp, 0, np.pi/20, '0', describer=None)
        simVertq, simVertI, simVerterr = ansect.sector(
            sim, 0, np.pi/20, '0', describer=None)

        expHorizq, expHorizI, expHorizerr = ansect.sector(
            exp, np.pi/2, np.pi/20, '0', describer=None)
        simHorizq, simHorizI, simHorizerr = ansect.sector(
            sim, np.pi/2, np.pi/20, '0', describer=None)

        expq, expI, experr = ansect.sector(
            exp, np.pi/2, np.pi, '0', describer=None)
        simq, simI, simerr = ansect.sector(
            sim, np.pi/2, np.pi, '0', describer=None)

        if pltErr == 1:
            fig = plt.figure(figsize=[8.5, 3], dpi=150)
            ax1 = fig.add_subplot(1, 3, 1)
            ax1.errorbar(expVertq, expVertI, yerr=expVerterr,
                         marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                         label='dataVert')
            ax1.plot(simVertq, simVertI,
                     marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                     label='simVert')
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            # print(expVerterr)
            ax1.legend()
            ax1.minorticks_off()

            ax2 = fig.add_subplot(1, 3, 2)
            ax2.errorbar(expHorizq, expHorizI, yerr=expHorizerr,
                         marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                         label='dataHoriz')
            ax2.plot(simHorizq, simHorizI, marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='',
                     label='simHoriz')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.legend()
            ax2.minorticks_off()

            ax3 = fig.add_subplot(1, 3, 3)
            ax3.errorbar(expq, expI, yerr=experr,
                         marker='o', markersize=2, markerfacecolor=[0, 0, 0], linestyle='', label='dataRad')
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

            ax2 = fig.add_subplot(1, 3, 2)
            ax2.plot(expHorizq, expHorizI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='dataHoriz')
            ax2.plot(simHorizq, simHorizI, marker='', markersize=2, markerfacecolor=[0, 0, 0], linestyle='-',
                     label='simHoriz')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.legend()
            ax2.minorticks_off()

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

        return

    def annularCompPlot(self, exp, sim, pltErr=0):
        """Plot annulus of experimental and modelled data.

        Args:
            exp: 2D sasmodels object of experimental data
            sim: 2D sasmodels object of model data

        Returns:
            Plot of experiment and model annulus
        """

        expVertq, expVertI, expVerterr = ansect.annular(exp, radius=0.07, thx=0.01)
        simVertq, simVertI, simVerterr = ansect.annular(sim, radius=0.07, thx=0.01)

        fig = plt.figure(figsize=[3.64*2, 2.48*2])
        ax = fig.add_axes([1, 1, 1, 1])

        ax.scatter(expVertq, expVertI, marker='o', label='data')
        ax.scatter(simVertq, simVertI, marker='o', label='sim')
        ax.legend()
        ax.minorticks_off()

        plt.rc('axes', linewidth=1.5)
        plt.rc('axes', grid=False)
        plt.rc('axes', labelsize='small')

        # Font
        plt.rc('font', family='sans-serif')
        plt.rc('font', weight='normal')
        plt.rc('font', size=14)

    def surfacePlot(self, zzq):
        """3D plot.

        Args:
            zzq: grid of interpolated intensity

        Returns:
            surface plot
        """
        plt.rcParams["figure.figsize"] = 12.8, 9.6
        # Normalize the colors based on Z value
        norm = plt.Normalize(zzq.min(), zzq.max())
        colors = cm.jet(norm(zzq))
        ax = plt.axes(projection='3d')
        surf = ax.plot_surface(self.xxq, self.yyq, zzq, facecolors=colors, shade=False)
        surf.set_facecolor((0, 0, 0, 0))

        return

    def objFunc(self, p_guess, p_list):
        """Objective function minimised during optimisation.

        This is defined as the reduced chi**2 statistic between experimental data and
        modelled data. Uses Levenberg--Marquardt algorithm by default as defined in code.

        Args:
            p_guess: List of initial guesses
            p_lsit: list of fitting parameters

        Returns:
            reduced chi**2 array
        """
        # =============================================================================

        dof = len(self.expData.q_data) - len(p_guess)

        # List of parameters that must be consistent between both bands
        omni_pars = ['scale', 'background', 'radius', 'radius_pd',
                     'radius_effective', 'volfraction']

        if 'bandVal' not in p_list and self.bandVal == 1:  # don't need static
            self.pars.update(dict(zip(p_list, p_guess)))
            sim = self.calculator(**self.pars)

            print(dict(zip(p_list, p_guess)))

        elif any(item in omni_pars for item in p_list):  # need to update static
            # fitting over bandval, need to remove from dicitonary before updating pars
            if 'bandVal' in p_list:
                self.pars.update(dict(zip(p_list, p_guess)))
                del self.pars['bandVal']

                hold = {items: p_guess[idx] for idx, items in enumerate(p_list)}
                staticDict = {}
                for keys, vals in hold.items():
                    if keys in omni_pars:
                        staticDict.update({keys: vals})
                    else:
                        continue

                self.staticPars.update(staticDict)

                sim = p_guess[p_list.index('bandVal')]*self.calculator(**self.pars) + \
                    (1 - p_guess[p_list.index('bandVal')])*self.calculator(**self.staticPars)

                print(dict(zip(p_list, p_guess)))

            # not fitting over bandval, use object defined value
            elif 'bandVal' not in p_list:
                #                print('option 3')
                self.pars.update(dict(zip(p_list, p_guess)))

                hold = {items: p_guess[idx] for idx, items in enumerate(p_list)}
                staticDict = {}
                for keys, vals in hold.items():
                    if keys in omni_pars:
                        staticDict.update({keys: vals})
                    else:
                        continue

                self.staticPars.update(staticDict)

                sim = self.bandVal*self.calculator(**self.pars) + \
                    (1 - self.bandVal)*self.calculator(**self.staticPars)

        else:  # only fitting pars present, no need to update static
            #            print('option 4')
            try:
                self.simStatic
            except:
                self.simStatic = self.calculator(**self.staticPars)

            if 'bandVal' in p_list:

                self.pars.update(dict(zip(p_list, p_guess)))
#                simShear = self.calculator(**self.pars)

                sim = p_guess[p_list.index('bandVal')]*self.calculator(**self.pars) + \
                    (1 - p_guess[p_list.index('bandVal')])*self.simStatic
                # print(dict(zip(p_list, p_guess)))

            elif 'bandVal' not in p_list:

                self.pars.update(dict(zip(p_list, p_guess)))
#                simShear = self.calculator(**self.pars)

                sim = self.bandVal*self.calculator(**self.pars) + \
                    (1 - self.bandVal)*self.simStatic

        err = ((((self.expData.data - sim)**2 / self.expData.err_data**2)))/dof
#        print(np.nansum((  ((  (self.expData.data - sim)**2 / self.expData.err_data**2)))))
        print(np.nansum(err))

        return err

# =============================================================================
# End defintion of classes
# =============================================================================

#############################################################################


def input_sample_check(conc, shear, id):
    """Check that conc and shear match with that defined by supplied index"""
    path = '../'
    file = '/expDataFiles.csv'

    os.path.isfile(path + file)
    exp = pd.read_csv(path+file, index_col='index')

    shear_raw = exp.loc[id, 'shear']
    sample_raw = exp.loc[id, 'sample']

    for idx, char in enumerate(shear_raw):
        if char != ' ':
            continue
        else:
            shearIdx = idx
    #            print(char)
            break

    for idx, char in enumerate(sample_raw):
        if char != 'w':
            continue
        else:
            sampleIdx = idx
    #            print(char)
            break

    # print(shearIdx)
    x = shear_raw[0:shearIdx]
    y = sample_raw[0:sampleIdx]

    if conc != y or shear != x:
        raise Exception('Conc or shear mismatch')
    else:
        print('Sample input matches')

    return

#############################################################################


def build_save_location(conc, shear_f):
    """Build path string to folder where model data and statistics should be saved."""
    if shear_f == '0':
        location = '../2D_simFits/ReducedChi2_fits/{conc}wt_CAPB_SLES_2wt_NaCl/{conc}wt_static/'.format(
            conc=conc)
    else:
        location = '../2D_simFits/ReducedChi2_fits/{conc}wt_CAPB_SLES_2wt_NaCl/{conc}wt_{shear}ps/'.format(
            conc=conc, shear=shear_f)

    return location

#############################################################################


def fitInput(fitChoose):
    """Build p_list from fitChoose"""
    p_list = []
    for keys, values in fitChoose.items():
        if values == 1:
            p_list.append(keys)
        else:
            continue
    return p_list

#############################################################################


def files_list_reduce(filename, fieldnames):
    """ Create array of input reduction settings """

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


def extract_sector(sans, arg1, arg2, description):
    """Extract sector for experimental and modelled data and calculate chi_2.

    Args:
        sans: 2D sasmodels data object
        arg1: centre of sector
        arg2: width of sector
        description: Describe the sector, eg. "vertical"

    Returns:
        string containing chi_2 and description
    """
    q_exp, I_exp, err_exp = ansect.sector(sans.expData, arg1, arg2)
    q_sim, I_sim, err_sim = ansect.sector(sans.simData, arg1, arg2)

    chi_2 = str(round(np.nansum((((I_exp-I_sim)/err_exp)**2)/len(q_exp)), sans.dp))
    stat = 'reduced chi2 of ' + description + ' sector: ' + chi_2

    return stat

#############################################################################


def extract_annulus(sans, arg1, arg2, description):
    """Extract annulus for experimental and modelled data and calculate chi_2.

    Args:
        sans: 2D sasmodels data object
        arg1: radius of annulus
        arg2: thickness of annulus
        description: Describe the annulus

    Returns:
        string containing chi_2 and description
    """

    q_exp, I_exp, err_exp = ansect.annular(sans.expData, radius=arg1, thx=arg2)
    q_sim, I_sim, err_sim = ansect.annular(sans.simData, radius=arg1, thx=arg2)

    chi_2 = str(round(np.nansum((((I_exp-I_sim)/err_exp)**2)/len(q_exp)), sans.dp))
    stat = 'reduced chi2 of ' + description + ': ' + chi_2

    return stat

#############################################################################


def save(sans, describer, minParams, minPars, stats, location, fitInfo, description):
    """Save modelled data and statistics

    Args:
        sans: 2D sasmodels object
        describer: String to be added to file name
        minParams: Full parameter dictionary
        minPars: Values of fitted parameters
        stats: Fitting statistics
        location: File path where files will be saved
        fitInfo: Information about termination condition and fitting algorithm used

    Returns:
        .dat file of modelled scattering data in qx, qy, I(q) format
        .txt file of input parameters and statistics corresponding to minimum
    """

    while path.exists(location) == False:
        print('error: file path does not exist. Please input a valid file path')
        location = input('file path: ')

    for idx, char in enumerate(sans.expData.shear[0]):
        if char != ' ':
            continue
        else:
            shearIdx = idx
            break

    # Build name for modelled scattering data
    shear = sans.expData.shear[0][0:shearIdx]

    name = sans.expData.sample[0] + '_' + shear + 'ps'
    post1 = '_sim'
    type1 = '.dat'

    saveName1 = name + post1 + describer + '_'
    # versionNum1 = input("Input a version number: ")
    versionNum1 = description

    # Write modelled scattering data to 3 column dat file
    write_3_column(location + saveName1 + versionNum1 + type1, sans)

    # Build name for modelled scattering data statistics
    post2 = '_simInfo'
    type2 = '.txt'

    saveName2 = name + post2 + describer + '_'

    output = []

    # Build output file
    output.append('qmin = ' + str(sans.qmin))
    output.append('ftol = ' + str(fitInfo[0]))
    output.append('method = ' + str(fitInfo[1]))
    output.append(' ')

    for key, val in minParams.items():
        if type(val) == str:
            output.append(str(key) + '=' + str(val) + ',')
        else:
            output.append(str(key) + '=' + str(round(val, sans.dp)) + ',')
    output.append(' ')

    output.append(' static parameters ')
    for key, val in sans.staticPars.items():
        if type(val) == str:
            output.append(str(key) + '=' + str(val) + ',')
        else:
            output.append(str(key) + '=' + str(round(val, sans.dp)) + ',')

    output.append(' ')

    output.append('Fitting_performed_over_the_following_parameters:')
    for key in minPars.keys():
        output.append(str(key))

    output.append('Returned_the_following_goodness_of_fit_measures:')
    output = output + stats
    output.append(str(datetime.datetime.now()))

    # Write output to txt file
    with open(location + saveName2 + versionNum1 + type2, 'w') as file:
        for lines in output:
            file.write(lines)
            file.write("\n")

    print('file was saved with filename: ' + saveName1 + versionNum1 + type1)
    return

#############################################################################


def write_3_column(filename, sans):
    """Write three column txt file.

    Args:
        filename: full file path to where the file will be saved
        sans: 2D sasmodels data object

    Returns:
        .dat file of modelled scattering data in qx, qy, I(q) format
    """
    with open(filename, 'wt') as fh:
        fh.write("Qx  Qy  I\n")
        for x, y, z in zip(sans.expData.qx_data, sans.expData.qy_data, sans.optimSim):
            fh.write("%g  %g  %g\n" % (x, y, z))

#############################################################################


def loadDict(build, popPars):
    """Load parameter dictionaries from saved txt file.

    **May break depending on when the txt file was saved due to formatting differences.

    Args:
        build: File path to txt file
        popPars: List of parameters to remove

    Returns:
        pars and pars_static dictionaries
    """
    data = pd.read_csv(build,
                       sep='=', skiprows=4, nrows=27, header=None)
    pars = {}
    vars = []
    vals = []
    for id in data.iterrows():
        if 'type' in id[1][0]:
            vars.append(id[1][0])
            vals.append(id[1][1][0:-1])
        else:
            vars.append(id[1][0])
            vals.append(float(id[1][1][0:-1]))

    pars = (dict(zip(vars, vals)))
    data_static = pd.read_csv(build,
                              sep='=', skiprows=33, nrows=26, header=None)

    pars_static = {}
    vars = []
    vals = []
    for id in data_static.iterrows():
        if 'type' in id[1][0]:
            vars.append(id[1][0])
            vals.append(id[1][1][0:-1])
        else:
            vars.append(id[1][0])
            vals.append(float(id[1][1][0:-1]))

    pars_static = (dict(zip(vars, vals)))
    print(popPars)

    for par in popPars:
        pars.pop(par)
        pars_static.pop(par)

    return pars, pars_static
    #############################################################################
