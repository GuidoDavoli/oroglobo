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
import utils.filtering_ecmwf as orofilt_ecmwf
import oroglobo_plotting as oroplot
from scipy import signal
from scipy import ndimage

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

d=int(oropar.filtering_Ndx_ecmwf["d"])
Ndx=int(oropar.filtering_Ndx_ecmwf["Ndx"])

# Open the file with the model mean orography 
data_model_grid_orog = xr.open_dataset(path_data_in+netcdf_model_grid_orog_in)
lat_model = data_model_grid_orog.latitude.values
lon_model = data_model_grid_orog.longitude.values
nlat=len(lat_model)
nlon=len(lon_model)

############ FILTERING 


#### FILTERING ALONG axis 0 (at each lon, filter along MERIDIONAL DIRECTION)

orogsmooth0=np.zeros((nlat,nlon))

for j in range(nlon):
    
    D=Ndx
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    orogsmooth0[:,j]=ndimage.convolve1d(data_model_grid_orog.elev[:,j],F1,mode='mirror')  # the edges are the poles --> approximate with "mirror"

print('meridional filtering: done!')

#### FILTER THE PREVIOUS FIELD ALONG AXIS 1 (at each LAT, filter along ZONAL DIRECTION)

operational_orog_on_model_grid=np.zeros((nlat,nlon))

for i in range(nlat):
    
    D=Ndx
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    operational_orog_on_model_grid[i,:]=ndimage.convolve1d(orogsmooth0[i,:],F1,mode='wrap')  # the edges are datelines or greenwich --> wrap around the globe
    
print('zonal filtering: done!')


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

