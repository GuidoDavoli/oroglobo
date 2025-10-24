#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

this file contains routines for filtering orography 
following ECMWF methodology

"""

import numpy as np
import math



def smoother_ECMWF(r, d, D):
    
    R = np.absolute(r)
    
    a = D/2 - d
    b = D/2 + d
    
    if R < a:
        
        h = 1/D
        
    if a <= R and R <= b:
        
        h = 1/(2*D) + 1/(2*D)*np.cos( np.pi*(r-D/2+d) )/(2*d)
        
    if R > b:
        
        h = 0
        
    return h


def distance(a,b):
    
    dist=np.sqrt(a**2 + b**2)
    
    return dist


def filter_ECMWF_2D(d,D):


    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
    #### prepare the 2d square array representing the filter (scale in pixel)
    #### the filter is a square with an odd number of pixels (L)
    #### but L is determined in different ways depending on D (even or odd)
    
    if (D%2)==0: # D is even
    
        L=int((D/2+d)*2+1)
        
    else:
        
        L=int(math.ceil(D/2+d)*2+1) # ceil(x) = round to the smallest integer higher or equal to x 
        
    
    filt = np.zeros((L,L))
    
    i0_filt=-int((L-1)/2)
    j0_filt=-int((L-1)/2)
    
    for i in range(L):
        for j in range(L):
            
            dist=distance(i0_filt+i, j0_filt+j)
            filt[i,j]=smoother_ECMWF(dist, d, D)
            
    filt=filt/np.sum(filt)  #### normalization after filtering
                            #### see https://medium.com/@bdhuma/6-basic-things-to-know-about-convolution-daef5e1bc411
                            #### the idea is to put a normalization factor in front of the filter matrix
    
    return filt


def filter_ECMWF_1D(d,D):


    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
    #### prepare the 1d array representing the filter (scale in pixel)
    #### the filter has an odd number of pixels (L)
    #### but L is determined in different ways depending on D (even or odd)
    
    if (D%2)==0: # D is even
    
        L=int((D/2+d)*2+1)
        
    else:
        
        L=int(math.ceil(D/2+d)*2+1) # ceil(x) = round to the smallest integer higher or equal to x 
        
    
    filt = np.zeros(L)
    
    i0_filt=-int((L-1)/2)
    
    for i in range(L):
            
        dist=i0_filt+i
        filt[i]=smoother_ECMWF(dist, d, D)
            
    filt=filt/np.sum(filt)  #### normalization after filtering
                            #### see https://medium.com/@bdhuma/6-basic-things-to-know-about-convolution-daef5e1bc411
                            #### the idea is to put a normalization factor in front of the filter matrix
    
    return filt