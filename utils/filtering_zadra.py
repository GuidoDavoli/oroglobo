#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:19:46 2024

@author: guidodavoli

THIS MODULE CONTAINS ROUTINES NEEDED TO PERFORM OROGRAPHIC FILTERING
AS EXPLAINED IN ZADRA, 2018: 
    Notes on the new low-pass filter for the orography field
    https://collaboration.cmc.ec.gc.ca/science/rpn/drag_project/documents/topo_lowpass_filter.pdf

WHEN POSSIBLE, THE EXECUTION OF SINGLE FUNCTIONS IS OPTIMIZED USING NUMBA

"""

import numpy as np
from numba import jit
import scipy.ndimage as ndimage
import sys


@jit(nopython=True, parallel=True) # Set "nopython" mode for best performance, equivalent to @njit
def LowPassFilter_c1(rc):
    
    rc=int(rc)
    
    c1=2/rc
    
    return c1


@jit(nopython=True, parallel=True) # Set "nopython" mode for best performance, equivalent to @njit
def LowPassFilter_coeffs(m,rc,p):
    
    # this function returns c(m) as defined in eq 30. in the paper the name of the coefficient is "m+1" but is just a name, c is a function of m
    # m must be between 1 and p-1, as in the paper eq. 28 (Summation part)
    # m, rc, p must be integers ---> they are forced to int
    
    m=int(m)
    rc=int(rc)
    p=int(p)
    
    c=-999
    
    if m<1 or m>(p-1):
        
        print("WARNING: WRONG M VALUED PASSED TO THE LP FILTER COEFFICENTS FUNCTION")
        
    if m>=1 and m<=(p-1):
        
        a=np.sin(2*np.pi*m/rc)
        b=2*np.pi*m/rc
        c=np.sin(2*np.pi*m/p)
        d=2*np.pi*m/p
        
        c=(2/rc)*(a/b)*(c/d) 
    
    return c
        

@jit(nopython=True, parallel=True) # Set "nopython" mode for best performance, equivalent to @njit
def LowPassFilter_Sp(rc,p):
    
    # this function returns Sp as defined in eq 32.
    # m must be between 1 and p-1, as in the paper eq. 28 (Summation part)
    # rc, p must be integers ---> they are forced to int

    rc=int(rc)
    p=int(p)

    summation=0
    
    for m in range(1,p): # range(1,p) gives numbers from 1 to p-1

        summation+=LowPassFilter_coeffs(m, rc, p)
        
    Sp=LowPassFilter_c1(rc) + 2*summation

    return Sp


def LowPassFilter_axis0_UPDATE(data_2D,rc,p):
    
    #### AGGIORNA DOC
    # this is the 1D filter of eq. 28
    # with the normalized coefficients "c hat" (eq. 31)
    # m and n are not taken as input because they represent the positions of the "image pixels";
    # the filter automatically scan all the 1d array and returns the filtered array.
    # rc, p must be integers ---> they are forced to int
    #
    #
    # rc: is the cut-off scale parameter, which indicates the threshold wavelength 
    #     (as a number of "data_1D" pixels) beyond which the amplitude of the signal should be reduced 
    #     see figure 1 in the paper.
    #
    # p : is the truncation parameter, which basically controls the sharpness 
    #     of the filter: the larger the value of p, the sharper the filter 
    #     see figure 2 in the paper. IN THIS VERSION IS AN ARRAY OF LENGTH EQUAL
    #     TO DATA2D_AXIS0; EACH VALUE IS THE VALUE OF THE TRUNC PARAM AT EACH "ROW" OF THE ARRAY
    #
    #     IMPORTANT:
    #     practically, p is the half-length of the averaging window:
    #     the value of the "central" pixel will be combined with those of (p − 1)
    #     “neighbors to the right” and (p − 1) “neighbors to the left”
    #
    #     THEREFORE,
    #     this function DO NOT FILTERS ALL PIXELS IN THE ARRAY; it excludes
    #     the firsts and lasts (p-1) pixels.
    #     The function returns an array with the shape of data_1D, with the
    #     "internal" pixel effectively filtered and the external pixels 
    #     (not filtered) ALL SET TO ZERO.
    
    rc=rc.astype('int')
    p=p.astype('int')
    
    datalen0=data_2D.shape[0]
    
    if len(p)!=datalen0:
        print("ERROR: len(p)!=datalen0",len(p),datalen0)
        sys.exit()
        
    if len(rc)!=datalen0:
        print("ERROR: len(rc)!=datalen0",len(rc),datalen0)
        sys.exit()
    
    data2D_filtered0=np.zeros(data_2D.shape)
    
    pixelcounter=1    
    
    for n in range(datalen0):
        
        print("n",n)
        
        if rc[n]>0: # if filtering (not in the padding region where rc=0)
            
            filter_window_half_len=p[n]-1
            nofilt_left_border=filter_window_half_len
            nofilt_right_border=datalen0-filter_window_half_len
                    
            if pixelcounter>nofilt_left_border and pixelcounter<nofilt_right_border: # inside the filtering region
            
                summation=0
                
                for  m in range(1,p[n]): # range(1,p) gives numbers from 1 to p-1
                
                    summation+=LowPassFilter_coeffs(m, rc[n], p[n])/LowPassFilter_Sp(rc[n], p[n]) * (data_2D[n+m,:] + data_2D[n-m,:])
                
                data2D_filtered0[n,:]=LowPassFilter_c1(rc[n])*data_2D[n,:] + summation
            
        pixelcounter+=1
    
    return data2D_filtered0


def LowPassFilter_axis1_UPDATE(data_2D,rc,p):
    
    #### AGGIORNA DOC
    # this is the 1D filter of eq. 28
    # with the normalized coefficients "c hat" (eq. 31)
    # m and n are not taken as input because they represent the positions of the "image pixels";
    # the filter automatically scan all the 1d array and returns the filtered array.
    # rc, p must be integers ---> they are forced to int
    #
    #
    # rc: is the cut-off scale parameter, which indicates the threshold wavelength 
    #     (as a number of "data_1D" pixels) beyond which the amplitude of the signal should be reduced 
    #     see figure 1 in the paper.
    #
    # p : is the truncation parameter, which basically controls the sharpness 
    #     of the filter: the larger the value of p, the sharper the filter 
    #     see figure 2 in the paper. IN THIS VERSION IS AN ARRAY OF LENGTH EQUAL
    #     TO DATA2D_AXIS1; EACH VALUE IS THE VALUE OF THE TRUNC PARAM AT EACH "COLUMN" OF THE ARRAY
    #
    #     IMPORTANT:
    #     practically, p is the half-length of the averaging window:
    #     the value of the "central" pixel will be combined with those of (p − 1)
    #     “neighbors to the right” and (p − 1) “neighbors to the left”
    #
    #     THEREFORE,
    #     this function DO NOT FILTERS ALL PIXELS IN THE ARRAY; it excludes
    #     the firsts and lasts (p-1) pixels.
    #     The function returns an array with the shape of data_1D, with the
    #     "internal" pixel effectively filtered and the external pixels 
    #     (not filtered) ALL SET TO ZERO.
    
    rc=rc.astype('int')
    p=p.astype('int')
    
    datalen1=data_2D.shape[1]
    
    if len(p)!=datalen1:
        print("ERROR: len(p)!=datalen1",len(p),datalen1)
        exit()
        
    if len(rc)!=datalen1:
        print("ERROR: len(rc)!=datalen1",len(rc),datalen1)
        exit()
    
    data2D_filtered1=np.zeros(data_2D.shape)
    
    pixelcounter=1    
    
    for n in range(datalen1):
        
        print("n",n)
        
        if rc[n]>0: # if filtering (not in the padding region where rc=0)
        
            filter_window_half_len=p[n]-1
            nofilt_left_border=filter_window_half_len
            nofilt_right_border=datalen1-filter_window_half_len
            
            if pixelcounter>nofilt_left_border and pixelcounter<nofilt_right_border: # inside the filtering region
            
                summation=0
                
                for  m in range(1,p[n]): # range(1,p) gives numbers from 1 to p-1
                
                    summation+=LowPassFilter_coeffs(m, rc[n], p[n])/LowPassFilter_Sp(rc[n], p[n]) * (data_2D[:,n+m] + data_2D[:,n-m])
                
                data2D_filtered1[:,n]=LowPassFilter_c1(rc[n])*data_2D[:,n] + summation
                
        pixelcounter+=1
    
    return data2D_filtered1



def LowPassFilter_2D_v2_UPDATE(data2D,rc0,rc1,p0,p1):
    
    # p0 and p1 are array of p values, they must have dimensions:
    # p0 equal to data2D.shape[0]
    # p1 equal to data2D.shape[1]
    
    # rc0 and rc1 are array of rc values, they must have dimensions:
    # rc0 equal to data2D.shape[0]
    # rc1 equal to data2D.shape[1]
    
    rc0=rc0.astype('int')
    rc1=rc1.astype('int')
    p0=p0.astype('int')
    p1=p1.astype('int')
    
    #data2D_filtered0=np.zeros(data2D.shape)
    #data2D_filtered=np.zeros(data2D.shape)

    data2D_filtered0=LowPassFilter_axis0_UPDATE(data2D, rc0, p0)
    data2D_filtered =LowPassFilter_axis1_UPDATE(data2D_filtered0, rc1, p1)
    
    return data2D_filtered


def localminimum_2D(data2D,N):
    
    # RETURN A MAP OF THE MINIMUM VALUE AMONG
    # THE NxN NEIGHBOURING POINTS (SEE ZADRA 2018, SECTION 5.2)
    #
    # CAUTION: NEED TO BETTER DEFINE TREATMENT OF BORDERS
    
    return ndimage.minimum_filter(data2D,N)


def localmaximum_2D(data2D,N):
    
    # RETURN A MAP OF THE MAXIMUM VALUE AMONG
    # THE NxN NEIGHBOURING POINTS (SEE ZADRA 2018, SECTION 5.2)
    #
    # CAUTION: NEED TO BETTER DEFINE TREATMENT OF BORDERS
    
    return ndimage.maximum_filter(data2D,N)


def apply_local_minmax_constraint(data2D_filtered,data2D_original,N):
    
    # APPLIES THE LOCAL MINIMUM/MAXIMUM CONSTRAINT AS DESCRIBED IN
    # ZADRA, 2018, SECTION 5.2. THE SIZE OF THE SQUARE DEFINING NEIGHBORING
    # POINTS IS PASSED BY THE USER (N)
    
    Hlmin=localminimum_2D(data2D_original,N)
    Hlmax=localmaximum_2D(data2D_original,N)
    
    return np.minimum( np.maximum( data2D_filtered,Hlmin ) , Hlmax )


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