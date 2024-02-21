#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:58:24 2023

@author: guidodavoli
"""

"""
THIS IS THE OROGLOBO PARAMETERS FILE / NAMELIST
in this files are defined variables used from other sections of the code 
parameters and parameter groups are defined through python dictionaries 

"""


paths_in = {
	"srtm30_data":"data_in/srtm30/",
	"model_grids":"data_in/model_grids/",
	"model_masks":"data_in/model_masks/",
    "copernicus_90m":"/home/guidodavoli/cnr_work/dati/Copernicus_DEM/copernicus_dem_90m_*SECTOR*/"
}

paths_out = {
	"img_out":"img_out/*GRIDNAME*/",
    "data_out":"data_out/*GRIDNAME*/"
}

paths_work = {
	"workdir":"work/",
    "workdir_grid":"work/*GRIDNAME*/"
}

files_in = {
	"netcdf_srtm30_global_nan_to_zero_0_360":"srtm30_global_nan_to_zero_0_360.nc",
	"netcdf_model_mask":"masks_*GRIDNAME*.nc",
	"txt_model_grid_lat":"*GRIDNAME*_lat.txt",
	"txt_model_grid_lon":"*GRIDNAME*_lon.txt",
    "tif_copernicus_90m":"Copernicus_DSM_30_*HEMISPHERE-NS**LAT2DIGITS*_00_*HEMISPHERE-EW**LON3DIGITS*_00_DEM.tif"	
}

files_out = {
	"img_1km_global_orog":"1km_global_orog.png",
	"img_1km_global_orog_smooth":"1km_global_orog_smooth.png",
    "img_model_grid_orog":"model_grid_orog_*GRIDNAME*.png",
    "img_model_grid_operational_orog":"model_grid_operational_orog_*GRIDNAME*.png",
    "netcdf_model_grid_operational_orog":"model_grid_operational_orog_*GRIDNAME*.nc",
    "netcdf_model_grid_ogwd_params":"model_grid_ogwd_params_*GRIDNAME*.nc",
    "netcdf_model_grid_tofd_params":"model_grid_tofd_params_*GRIDNAME*.nc"
}

files_work = {
	"netcdf_1km_smooth_orog":"1km_global_orog_smooth_to_*GRIDNAME*.nc",
	"netcdf_model_grid_orog":"model_grid_orog_*GRIDNAME*.nc"
}

model_grid = {
	"GRIDNAME":"globo_KM312",
    "gridspacing_max_km":312
}




