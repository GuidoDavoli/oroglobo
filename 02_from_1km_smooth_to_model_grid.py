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

# START FROM THE FIRST LATITUDE --> SOUTH POLE

# GLOBO POINTS ARE THE T POINTS, SO POINTS AT LAT +90 AND -90 REPRESENT
# CIRCULAR AREAS AROUND THE POLES

# THE CODE EXPECT THE DATA ORDERED IN LATITUDE FROM -90 TO +90

mean_orog_on_model_grid=np.zeros((len(lat_model),len(lon_model)))

for i in range(len(lat_model)):
    
    if lat_model[i]==-90:
        
        print("SP")
        
        dr=abs(lat_model[i]-lat_model[i+1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        northboundary=-90+dr # lat boundary of the polar cap
        
        # collect all the "1km points" falling in the polar cap
        
        polar_cap_points_1km=data_orog_1km_smooth.elev.loc[data_orog_1km_smooth.latitude<northboundary]
        
        mean_orog_in_this_model_gridpoint=polar_cap_points_1km.data.mean()
        
        mean_orog_on_model_grid[i,:]=mean_orog_in_this_model_gridpoint # all points at -90 have the same value, at every longitude
        
    if lat_model[i]==90:
        
        print("NP")
        
        dr=abs(lat_model[i]-lat_model[i-1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        southboundary=90-dr # lat boundary of the polar cap
        
        # collect all the "1km points" falling in the polar cap
        
        polar_cap_points_1km=data_orog_1km_smooth.elev.loc[data_orog_1km_smooth.latitude>southboundary]
        
        mean_orog_in_this_model_gridpoint=polar_cap_points_1km.data.mean()
        
        mean_orog_on_model_grid[i,:]=mean_orog_in_this_model_gridpoint # all points at +90 have the same value, at every longitude

        
    if lat_model[i]>-90 and lat_model[i]<90:
        
        print(lat_model[i])
        
        for j in range(len(lon_model)):
            
            print(lon_model[j])
