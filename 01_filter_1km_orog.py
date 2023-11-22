#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:58:16 2023

@author: guido
"""

import numpy as np
import xarray as xr
import oroglobo_parameters as oropar
import oroglobo_functions as orofunc
import oroglobo_plotting as oroplot
from scipy import signal


# IMPORT PARAMETERS
gridname=oropar.model_grid["GRIDNAME"]

path_data_in=oropar.paths_in["srtm30_data"]
path_img_out=oropar.paths_out["img_out"].replace("*GRIDNAME*", gridname) 
path_data_out=oropar.paths_work["workdir_grid"].replace("*GRIDNAME*", gridname) 

netcdf_orog_in=oropar.files_in["netcdf_srtm30_global_nan_to_zero_0_360"]

img_1km_global_raw_out=oropar.files_out["img_1km_global_orog"]
img_1km_global_smooth_out=oropar.files_out["img_1km_global_orog_smooth"]
netcdf_1km_smooth_orog_out=oropar.files_work["netcdf_1km_smooth_orog"].replace("*GRIDNAME*", gridname) 



# Open a file
data_orog = xr.open_dataset(path_data_in+netcdf_orog_in)
lat = data_orog.latitude.values
lon = data_orog.longitude.values


# plot
oroplot.orography_plot(data_orog.elev,path_img_out+img_1km_global_raw_out,2400)


############ FILTERING 


filt=orofunc.filter_ECMWF(1, 80)

orogsmooth = signal.convolve2d(data_orog.elev, filt, mode='same', boundary='symm').astype(np.float32) # IF BUNDARY IS SET TO WRAP, SPURIOUS OROGRAPHY AT THE POLES (np IS WRAPPED WIT sp AND VICE VERSA)
# ANOTHER SOLUTION WOULD BE TO USE "WRAP" TO GET CORRECT WRAPPING ALONG LONGITUDES, ADN CORRECT AT NP and SP: 0 at NP points, a simple average at SP points

# plot
oroplot.orography_plot(orogsmooth,path_img_out+img_1km_global_smooth_out,2400)

#orogsmooth_da=xr.DataArray(orogsmooth.transpose, coords=data_orog.coords, dims=data_orog.dims, attrs=data_orog.attrs)
orogsmooth_da=xr.DataArray(orogsmooth, coords=[('latitude', lat),('longitude', lon)])

# save as netcdf
orogsmooth_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_smooth_orog_out)


# su tintin gira in 4 ore

## WARNING: AT PRESENT THE SMOOTHING OPERATION CREATE SPURIOUS OROGRAPHY OVER NORTH POLE