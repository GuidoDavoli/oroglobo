#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 10:19:21 2023

@author: guidodavoli

THIS SCRIPT CREATE THE OPERATIONAL GRID:
    
    1) SMOOTHS THE MODEL GRIDBOX-AVERAGED OROGRAPHY
    2) APPLIES THE MODEL LAND-SEA MASK (ASSUMED LON ORDER 0-360)

"""

import numpy as np
import xarray as xr
import oroglobo_parameters as oropar
import oroglobo_functions as orofunc
import oroglobo_plotting as oroplot
from scipy import signal

# IMPORT PARAMETERS
gridname=oropar.model_grid["GRIDNAME"]

path_data_in=oropar.paths_work["workdir_grid"].replace("*GRIDNAME*", gridname) 
path_data_in_mask=oropar.paths_in["model_masks"].replace("*GRIDNAME*", gridname) 
path_data_out=oropar.paths_out["data_out"].replace("*GRIDNAME*", gridname) 
path_img_out=oropar.paths_out["img_out"].replace("*GRIDNAME*", gridname) 

netcdf_model_grid_orog_in=oropar.files_work["netcdf_model_grid_orog"].replace("*GRIDNAME*", gridname) 
netcdf_model_mask_in=oropar.files_in["netcdf_model_mask"].replace("*GRIDNAME*", gridname) 
netcdf_model_grid_operational_orog_out=oropar.files_out["netcdf_model_grid_operational_orog"].replace("*GRIDNAME*", gridname) 
img_model_grid_operational_orog_out=oropar.files_out["img_model_grid_operational_orog"].replace("*GRIDNAME*", gridname) 



# Open the file with the model mean orography 
data_model_grid_orog = xr.open_dataset(path_data_in+netcdf_model_grid_orog_in)
lat_model = data_model_grid_orog.latitude.values
lon_model = data_model_grid_orog.longitude.values

############ FILTERING 

filt=orofunc.filter_ECMWF(1, 6)

operational_orog_on_model_grid = signal.convolve2d(data_model_grid_orog.elev, filt, mode='same', boundary='symm').astype(np.float32) # IF BUNDARY IS SET TO WRAP, SPURIOUS OROGRAPHY AT THE POLES (np IS WRAPPED WIT sp AND VICE VERSA)
# ANOTHER SOLUTION WOULD BE TO USE "WRAP" TO GET CORRECT WRAPPING ALONG LONGITUDES, ADN CORRECT AT NP and SP: 0 at NP points, a simple average at SP points

############ MASKING

# Open the mask
data_model_mask = xr.open_dataset(path_data_in_mask+netcdf_model_mask_in)
land_sea_mask = data_model_mask["glog.msk"].values[:,:-2]

# Apply the mask
# THIS IS WRONG, PERCHÈ I LAGHI SONO =0 COME MASK MA IN REALTÀ NON SONO A ZERO COME ALTITUDINE
operational_orog_on_model_grid=np.where(land_sea_mask==1, operational_orog_on_model_grid, 0)
operational_orog_on_model_grid=np.where(land_sea_mask==0, 0, operational_orog_on_model_grid)


### PLOT
oroplot.orography_plot(operational_orog_on_model_grid,path_img_out+img_model_grid_operational_orog_out,600)


### create data array and save to netcdf

operational_orog_on_model_grid_da=xr.DataArray(operational_orog_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])

operational_orog_on_model_grid_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_model_grid_operational_orog_out)

