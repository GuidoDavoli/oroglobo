#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:58:16 2023

@author: guido
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import oroglobo_parameters as oropar

"""
SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
"""


def smoothing2D_ECMWF(r, d, D):
    
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


# IMPORT PARAMETERS
path_data_in=oropar.paths_in["srtm30_data"]
path_img_out=oropar.paths_out["img_out"]
path_data_out=oropar.paths_work["workdir"]

netcdf_orog_in=oropar.files_in["netcdf_srtm30_global_nan_to_zero_0_360"]

img_srtm30_global_nan_to_zero_out=oropar.files_out["img_srtm30_global_nan_to_zero"]
img_srtm30_smooth_out=oropar.files_out["img_srtm30_smooth"]
netcdf_1km_smooth_orog_out=oropar.files_work["netcdf_1km_smooth_orog"]



# Open a file
data_orog = xr.open_dataset(path_data_in+netcdf_orog_in)
lat = data_orog.latitude.values
lon = data_orog.longitude.values


# print
plt.figure()
plt.imshow(data_orog.elev,origin='lower',cmap='terrain',vmin=0,vmax=8500)
plt.colorbar(location='bottom')
plt.savefig(path_img_out+img_srtm30_global_nan_to_zero_out,dpi=2400)
plt.show()


############ FILTERING 

from scipy import signal

#### prepare the 2d square array representing the filter (scale in pixel)

d = 1
D = 80

L=2*int(D/2+2)

filt = np.zeros((L,L))

i0_filt=-L/2
j0_filt=-L/2

for i in range(L):
    for j in range(L):
        
        dist=distance(i0_filt+i, j0_filt+j)
        filt[i,j]=smoothing2D_ECMWF(dist, d, D)
        
filt=filt/np.sum(filt)  #### tentativo di conservare il "modulo" dell'orografia
                        #### ispirato da https://medium.com/@bdhuma/6-basic-things-to-know-about-convolution-daef5e1bc411
                        #### the idea is to put a normalization factor in front of the filter matrix

orogsmooth = signal.convolve2d(data_orog.elev, filt, mode='same', boundary='wrap').astype(np.float32)

plt.figure()
plt.imshow(orogsmooth,origin='lower',cmap='terrain',vmin=0,vmax=8500)
plt.colorbar(location='bottom')
plt.savefig(path_img_out+img_srtm30_smooth_out,dpi=2400)
plt.show()


#orogsmooth_da=xr.DataArray(orogsmooth.transpose, coords=data_orog.coords, dims=data_orog.dims, attrs=data_orog.attrs)
orogsmooth_da=xr.DataArray(orogsmooth, coords=[('latitude', lat),('longitude', lon)])

# save as netcdf
orogsmooth_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_smooth_orog_out)


# su tintin gira in 4 ore
