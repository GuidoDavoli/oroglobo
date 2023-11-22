#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:48:14 2023

@author: guidodavoli
"""

import numpy as np
import math


def distance(a,b):
    
    dist=np.sqrt(a**2 + b**2)
    
    return dist


def smoothing2D_ECMWF(r, d, D):
    
    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
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


def filter_ECMWF(d,D):


    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
    #### prepare the 2d square array representing the filter (scale in pixel)
    #### the filter is a square with an odd number of pixels (L)
    #### but L is determinig in different ways depending on D (even or odd)
    
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
            filt[i,j]=smoothing2D_ECMWF(dist, d, D)
            
    filt=filt/np.sum(filt)  #### tentativo di conservare il "modulo" dell'orografia
                            #### ispirato da https://medium.com/@bdhuma/6-basic-things-to-know-about-convolution-daef5e1bc411
                            #### the idea is to put a normalization factor in front of the filter matrix
    
    return filt


