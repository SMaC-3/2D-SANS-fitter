#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:04:38 2020

@author: jkin0004
"""

import numpy as np
import math
import rheoSANSFunctions_fitOpt_lsq as rsf
import matplotlib.pyplot as plt

def extractLoop():    
    indexSelect = '6,15,24,34,43,62,71,83,95'

    indexNums = rsf.evaluate_files_list(indexSelect)
    describer = 'radAve'
    saveSet = '1'
    
    centre = np.pi/2
    width = np.pi
    
    radius = 0.07
    thx = 0.01
    
    sans = rsf.sans2d()
    sans.qmin = 0.007
    
    for nums in indexNums:
        sans.getData(str(nums))
        
        q, I_q, I_err = sector(sans.expData, centre = centre, width = width, 
                                      save=saveSet, describer = describer)
        
#        q_ann, I_q_ann = annular(sans.expData, radius= radius,thx= thx,
#                                 save=saveSet, describer = describer)
        del sans.expData
        
        fig = plt.figure()
        
        kwargs = {'xscale' : 'log', 'yscale': 'log'}
        
        ax = fig.add_axes([1,1,1,1], **kwargs)
        
        ax.errorbar(q,I_q, yerr = I_err, linestyle = '', marker = 'o', markersize = 1)

    
    return

def sector(dataSet, centre, width, save='0', describer = None):
    
    xcentre=0
    ycentre=0
    
    sectcent= centre # In radians
    sectwid= width # In radians
    
    #dataSet should be a sasView 2D data object
    #data_sqrd = dataSet**2
    mag = dataSet.q_data #q values
    sectang = []

    
    for i in range(len(mag)):
        if dataSet.qy_data[i] > 0:
            sectang.append(math.atan(dataSet.qx_data[i]/dataSet.qy_data[i]))
        elif dataSet.qy_data[i] < 0:
            sectang.append(math.atan(dataSet.qx_data[i]/dataSet.qy_data[i]) + np.pi)
    
    sectang = np.array(sectang)
    
    cwmax=sectcent+(0.5*sectwid) #Max bound for top sector
    cwmin=sectcent-(0.5*sectwid) #Min bound for top sector
    cwmax_ref=sectcent+(0.5*sectwid)+math.pi #Max bound for bottom sector
    cwmin_ref=sectcent-(0.5*sectwid)+math.pi #Min bound for bottom sector
    
    #Sorting according to angle
    sortI = np.argsort(sectang)
    sectang = sectang[sortI]
    mag = mag[sortI]
    err = dataSet.err_data[sortI]
    data = dataSet.data[sortI]
    
    
    seccrop_data = np.zeros_like(data)
    seccrop_err = np.zeros_like(err)
    seccrop_mag = np.zeros_like(mag)
    seccrop_ang = np.zeros_like(sectang)
    
    #Find logic gates
    posLog = np.logical_and(cwmin<sectang, sectang<cwmax)
    negLog = np.logical_and(cwmin_ref<sectang, sectang<cwmax_ref)
    
    #Find values according to logic gates
    seccrop_data[posLog] = data[posLog]
    seccrop_err[posLog] = err[posLog]
    seccrop_mag[posLog] = mag[posLog]
    seccrop_ang[posLog] = sectang[posLog]
    
    seccrop_data[negLog] = data[negLog]
    seccrop_err[negLog] = err[negLog]
    seccrop_mag[negLog] = mag[negLog]
    seccrop_ang[negLog] = sectang[negLog]
    
    zeros = seccrop_mag != 0 #Find zeros
    
    #remove zeros
    seccrop_data = seccrop_data[zeros]
    seccrop_err = seccrop_err[zeros]
    seccrop_mag = seccrop_mag[zeros]
    seccrop_ang = seccrop_ang[zeros]
    
    #Sort by ascending q
    sortq = np.argsort(seccrop_mag)
    
    seccrop_data = seccrop_data[sortq]
    seccrop_err = seccrop_err[sortq]
    seccrop_mag = seccrop_mag[sortq]
    seccrop_ang = seccrop_ang[sortq]
    
    #Make 100 bins spaced log linearly between min and max q
    nbs = 100
    minMag = np.min(seccrop_mag)
    maxMag = np.max(seccrop_mag)
    
    logMinMag = np.log10(minMag)
    logMaxMag = np.log10(maxMag)
    
    logLinear = np.linspace(logMinMag, logMaxMag, nbs)
    bins = 10**logLinear
    
    binindex = np.zeros([nbs]) #number points summed per bin
    bintotal = np.zeros([nbs]) #summed intensity
    errtotal = np.zeros([nbs]) #summed error
    
    for i in range(nbs - 1):
        for ii in range(len(seccrop_data)):
            if seccrop_mag[ii]>=bins[i]:
                if seccrop_mag[ii] <= bins[i + 1]:
                    binindex[i] = binindex[i] + 1
                    bintotal[i] = bintotal[i] + seccrop_data[ii]
                    errtotal[i] = errtotal[i] + seccrop_err[ii]
#                print(errtotal[i])
                
    binZeros = binindex != 0 
    
    bins = bins[binZeros]
    binindex = binindex[binZeros]
    bintotal = bintotal[binZeros]
    errtotal = errtotal[binZeros]
#    print(errtotal)
    
    binave = bintotal/binindex
    errave = errtotal/binindex
    
#    allerror = [err, seccrop_err, errtotal, errave]
    
    if save == '1':
        fileType = '.dat'
        if dataSet.shear[0][0] == '0':
            fileName = describer + '_' + str(dataSet.sample[0]) + '_' + 'static'
        else:
            fileName = describer + '_' + str(dataSet.sample[0]) + '_' + str(dataSet.shear[0][0:-14]) + 'ps'
        location = '../2D_annular_sector_extraction/py_sect_radAve/'
        fullName = location + fileName + fileType
        with open(fullName, 'wt') as fh:
            fh.write("q  I(q)  err_I\n")
            for x, y, z in zip(bins, binave, errave):
                fh.write("%g  %g  %g\n" % (x, y, z))
    
    
    return bins, binave, errave

def annular(dataSet, radius, thx, save='0', describer = None):
    
    radius = radius
    thx = thx
    
    mag = dataSet.q_data
     
    #Draw a set of x,y points for the circles chosen
    theta = np.linspace(0, 2*np.pi, 314)
    cx = radius*np.cos(theta);
    cy = radius*np.sin(theta);
    cxouter = (radius + thx)*np.cos(theta);
    cyouter = (radius + thx)*np.sin(theta);
    
    #Capture points that fall within the annular rings based on their magnitude
    
    I_ann = np.logical_and(radius<mag , mag < (radius+thx))
    annul_x = dataSet.qx_data[I_ann]
    annul_y = dataSet.qy_data[I_ann]
    annul_I = dataSet.data[I_ann]
    annul_mag = dataSet.q_data[I_ann]
    
    #Calculate the angles for the obtained points. Zero is twleve o clock (vertically up on y-axis)

    annul_ang = []
    
    for i in range(len(annul_mag)):
        if annul_y[i] > 0:
            annul_ang.append(math.atan(annul_x[i]/annul_y[i]))
        elif annul_y[i] < 0:
            annul_ang.append(math.atan(annul_x[i]/annul_y[i]) + np.pi)
    
    annul_ang = np.array(annul_ang)
    
    #Sorts data to give as a function of increasing angle
    sortI = np.argsort(annul_ang)
    annul_ang = annul_ang[sortI]
    annul_mag = annul_mag[sortI]
    annul_x = annul_x[sortI]
    annul_y = annul_y[sortI]
    annul_I = annul_I[sortI]
    
    #Data binning
    nbsa = 100
    deltheta = 2*np.pi/nbsa
    binsa = np.linspace(-np.pi/2, 3*np.pi/2, nbsa)
    
    binindexa = np.zeros([nbsa]) #number points summed per bin
    bintotala = np.zeros([nbsa]) #summed intensity
#    errtotala = np.zeros([nbs]) #summed error
    
    for i in range(nbsa - 1):
        for ii in range(len(annul_mag)):
            if annul_ang[ii]>=binsa[i]:
                if annul_ang[ii] <= binsa[i + 1]:
                    binindexa[i] = binindexa[i] + 1
                    bintotala[i] = bintotala[i] + annul_I[ii]
#                    errtotal[i] = errtotal[i] + seccrop_err[ii]
#                print(errtotal[i])
                
    binZeros = binindexa != 0 
    
    binsa = binsa[binZeros]
    binindexa = binindexa[binZeros]
    bintotala = bintotala[binZeros]
#    errtotala = errtotala[binZeros]
#    print(errtotal)
    
    binavea = bintotala/binindexa
#    erravea = errtotala/binindexa
#    print(errave)
    
    if save == '1':
        fileType = '.dat'
        fileName = describer + '_' + str(dataSet.sample[0]) + '_' + str(dataSet.shear[0][0:-14]) + 'ps'
        location = '../2D_annular_sector_extraction/py_annular/'
        fullName = location + fileName + fileType
        with open(fullName, 'wt') as fh:
            fh.write("q  I(q)  err_I\n")
            for x, y, z in zip(binsa, binavea, binavea):
                fh.write("%g  %g  %g\n" % (x, y, z))
    
    return binsa, binavea