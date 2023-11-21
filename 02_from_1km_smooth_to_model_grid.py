#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:50:31 2023

@author: guidodavoli
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import oroglobo_parameters as oropar


# IMPORT PARAMETERS
path_grid_in=oropar.paths_in["model_grids"]
path_data_in=oropar.paths_work["workdir"]
path_img_out=oropar.paths_out["img_out"]
path_data_out=oropar.paths_work["workdir"]

txt_model_grid_lat_in=oropar.files_in["txt_model_grid_lat"]
txt_model_grid_lon_in=oropar.files_in["txt_model_grid_lon"]

netcdf_model_grid_orog_out=oropar.files_work["netcdf_model_grid_orog"]
netcdf_1km_smooth_orog_in=oropar.files_work["netcdf_1km_smooth_orog"]


# Open the file(s) with model gridpoints
lat_model=np.loadtxt(path_grid_in+txt_model_grid_lat_in)
lon_model=np.loadtxt(path_grid_in+txt_model_grid_lon_in)

# Open the file with 1km-resolution orography smoothed to model gridspacing
data_orog_1km_smooth = xr.open_dataset(path_data_in+netcdf_1km_smooth_orog_in)
lat_1km = data_orog_1km_smooth.latitude.values
lon_1km = data_orog_1km_smooth.longitude.values

# AVERAGE THE 1KM SMOOTH OROGRAPHY TO THE MODEL GRIDBOXES.
# AT PRESENT, THIS PART OF THE CODE IS 100% MODEL-DEPENDENT
# AND WORKS FOR GLOBO ONLY

check that model & 1km longitudes are 0 ... 360 and lats -90...90
 