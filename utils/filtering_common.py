#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:32:57 2024

@author: guidodavoli

COMMON ROUTINES FOR FILTERING IN OROGLOBO

"""

import scipy.ndimage as ndimage
import numpy as np
from multiprocessing import Pool


######### routines for local minmax constraint, for 2d data, with constant ####
######### neighborhood scale (in pixels) around the globe #####################

def localminimum_2D(data2D,N,how):
    
    # RETURN A MAP OF THE MINIMUM VALUE AMONG
    # THE NxN NEIGHBOURING POINTS (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    return ndimage.minimum_filter(data2D,N,mode=how)


def localmaximum_2D(data2D,N,how):
    
    # RETURN A MAP OF THE MAXIMUM VALUE AMONG
    # THE NxN NEIGHBOURING POINTS (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    return ndimage.maximum_filter(data2D,N,mode=how)


def apply_local_minmax_constraint_2D_constantN(data2D_filtered,data2D_original,N,how):
    
    # APPLIES THE LOCAL MINIMUM/MAXIMUM CONSTRAINT AS DESCRIBED IN
    # ZADRA, 2018, SECTION 5.2. THE SIZE OF THE SQUARE DEFINING NEIGHBORING
    # POINTS IS PASSED BY THE USER (N)
    # "how" specifies what to do with borders along the two axis
    
    if how[0]=='symmetric':
        how[0]='mirror'
        
    if how[1]=='symmetric':
        how[1]='mirror'
    
    Hlmin=localminimum_2D(data2D_original,N,how)
    Hlmax=localmaximum_2D(data2D_original,N,how)
    
    return np.minimum( np.maximum( data2D_filtered,Hlmin ) , Hlmax )

###############################################################################

######### routines for local minmax constraint, for 2d data, with variable ####
######### neighborhood scale (in pixels) around the globe #####################
######### (applied axis by axis) ##############################################


def localminimum_2D_axis1(data2D,Naxis1,how,Nparal=1):
    
    # RETURN A MAP OF THE local MINIMUM VALUE AMONG
    # axis1; lenght scale specified in Naxis1; can be variable (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    naxis0=data2D.shape[0]
    data2D_localmin_axis1=np.zeros(data2D.shape)
    
    for i in range(naxis0):
        
        data2D_localmin_axis1[i,:]=ndimage.minimum_filter1d(data2D[i,:],Naxis1[i],mode=how)
    """
    def funzA(i):
        data2D_localmin_axis1_i=ndimage.minimum_filter1d(data2D[i,:],Naxis1[i],mode=how) 
        return data2D_localmin_axis1_i
    
    pool = Pool(processes=Nparal)
    ####### HERE COMES THE CHANGE #######
    results = [pool.apply_async(funzA, [val]) for val in range(naxis0)]
    for idx, val in enumerate(results):
        data2D_localmin_axis1[idx,:]= val.get()
    #######
    pool.close()
    """
    return data2D_localmin_axis1


def localmaximum_2D_axis1(data2D,Naxis1,how,Nparal=1):
    
    # RETURN A MAP OF THE local MAXIMUM VALUE AMONG
    # axis1; lenght scale specified in Naxis1; can be variable (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    naxis0=data2D.shape[0]
    data2D_localmax_axis1=np.zeros(data2D.shape)
    
    for i in range(naxis0):
        
        data2D_localmax_axis1[i,:]=ndimage.maximum_filter1d(data2D[i,:],Naxis1[i],mode=how)
    """
    def funzB(i):
        data2D_localmax_axis1_i=ndimage.maximum_filter1d(data2D[i,:],Naxis1[i],mode=how)
        return data2D_localmax_axis1_i
    
    pool = Pool(processes=Nparal)
    ####### HERE COMES THE CHANGE #######
    results = [pool.apply_async(funzB, [val]) for val in range(naxis0)]
    for idx, val in enumerate(results):
        data2D_localmax_axis1[idx,:]= val.get()
    #######
    pool.close()
    """
    return data2D_localmax_axis1


def localminimum_2D_axis0(data2D,Naxis0,how,Nparal=1):
    
    # RETURN A MAP OF THE local MINIMUM VALUE AMONG
    # axis0; lenght scale specified in Naxis0; can be variable (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    naxis1=data2D.shape[1]
    data2D_localmin_axis0=np.zeros(data2D.shape)
    
    for j in range(naxis1):
        
        data2D_localmin_axis0[:,j]=ndimage.minimum_filter1d(data2D[:,j],Naxis0[j],mode=how)
    """
    def funzC(j):
        data2D_localmin_axis0_j=ndimage.minimum_filter1d(data2D[:,j],Naxis0[j],mode=how)
        return data2D_localmin_axis0_j
    
    pool = Pool(processes=Nparal)
    ####### HERE COMES THE CHANGE #######
    results = [pool.apply_async(funzC, [val]) for val in range(naxis1)]
    for idx, val in enumerate(results):
        data2D_localmin_axis0[:,idx]= val.get()
    #######
    pool.close()
    """
    return data2D_localmin_axis0


def localmaximum_2D_axis0(data2D,Naxis0,how,Nparal=1):
    
    # RETURN A MAP OF THE local MAXIMUM VALUE AMONG
    # axis0; lenght scale specified in Naxis0; can be variable (SEE ZADRA 2018, SECTION 5.2)
    # CALLED BY apply_local_minmax_constraint
    
    naxis1=data2D.shape[1]
    data2D_localmax_axis0=np.zeros(data2D.shape)
    
    for j in range(naxis1):
        
        data2D_localmax_axis0[:,j]=ndimage.maximum_filter1d(data2D[:,j],Naxis0[j],mode=how)
    """
    def funzD(j):
        data2D_localmax_axis0_j=ndimage.maximum_filter1d(data2D[:,j],Naxis0[j],mode=how)
        return data2D_localmax_axis0_j
    
    pool = Pool(processes=Nparal)
    ####### HERE COMES THE CHANGE #######
    results = [pool.apply_async(funzD, [val]) for val in range(naxis1)]
    for idx, val in enumerate(results):
        data2D_localmax_axis0[:,idx]= val.get()
    #######
    pool.close()
    """
    return data2D_localmax_axis0


def apply_local_minmax_constraint_2D_variableN(data2D_filtered,data2D_original,Nzonal,Nmeridional,how,Nparal=1):
    
    # APPLIES THE LOCAL MINIMUM/MAXIMUM CONSTRAINT AS DESCRIBED IN
    # ZADRA, 2018, SECTION 5.2. 
    # The constraint is applied axis by axis and N is variable
    # Naxis0 and Naxis1 have the dimensions of the two axis
    # "how" specifies what to do with borders along the two axis
    
    if how[0]=='symmetric':
        how[0]='mirror'
        
    if how[1]=='symmetric':
        how[1]='mirror'
    
    Hlmin0=localminimum_2D_axis0(data2D_original, Nmeridional, how[0], Nparal)
    Hlmin1=localminimum_2D_axis1(data2D_original, Nzonal, how[1], Nparal)
    Hlmin=np.minimum(Hlmin0,Hlmin1)
    del Hlmin0 # not needed anymore, free some memory
    del Hlmin1 # not needed anymore, free some memory

    Hlmax0=localmaximum_2D_axis0(data2D_original, Nmeridional, how[0], Nparal)
    Hlmax1=localmaximum_2D_axis1(data2D_original, Nzonal, how[1], Nparal)
    Hlmax=np.maximum(Hlmax0,Hlmax1)
    del Hlmax0 # not needed anymore, free some memory
    del Hlmax1 # not needed anymore, free some memory
    
    A=np.maximum( data2D_filtered,Hlmin )
    
    del data2D_filtered # not needed anymore, free some memory
    del Hlmin # not needed anymore, free some memory
    
    return np.minimum( A , Hlmax )




###############################################################################

######### generic routines for adding zero-padded frames around 2d data ####### 

def add_frame_before_filtering(data2D, D, how=['symmetric','symmetric']):
    
    # add a frame, with thickness D pixels, around data2D.
    # how to fill the frame along the 2 axis is decided with the "how" parameter
    # possible options for "how":
    #   - symmetric: the frame contains the same pixels as data2D, symmetrically wrt data2D border
    #   - wrap: the frame contains pixel from data2D taken from the opposite side of data2D
    # the 2 "how" options can be different for the 2 axis.
    
    # add pad along first axis
    
    temp = np.pad(data2D,((D,D),(0,0)),mode=how[0])
    
    # add pad along second axis
    
    data2D_with_frame = np.pad(temp,((0,0),(D,D)),mode=how[1])
    
    return data2D_with_frame


def remove_frame_after_filtering(data2D, D):
    
    # return an array without an external frame of thickness D
    
    ax0len=data2D.shape[0]
    ax1len=data2D.shape[1]
    
    return data2D[D:ax0len-D,D:ax1len-D]