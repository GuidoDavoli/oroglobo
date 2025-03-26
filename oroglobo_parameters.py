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
    "copernicus_90m":"/home/guidodavoli/cnr_work/dati/Copernicus_DEM/copernicus_dem_90m_*SECTOR*/",
    "copernicus_lowres":"/home/guidodavoli/cnr_work/dati/Copernicus_DEM/copernicus_dem_lowres/copernicus_dem_90m_*SECTOR*/",
    "copernicus_highresnc":"/home/guidodavoli/cnr_work/dati/Copernicus_DEM/copernicus_dem_highresnc/copernicus_dem_90m_*SECTOR*/"
}

paths_out = {
	"img_out":"img_out/*GRIDNAME*/",
    "data_out":"data_out/*GRIDNAME*/",
    "1km_out": "data_out/km1/"
}

paths_work = {
	"workdir":"work/",
    "workdir_grid":"work/*GRIDNAME*/"
}

files_in = {
	"netcdf_srtm30_global_nan_to_zero_0_360":"srtm30_global_nan_to_zero_0_360.nc",
	"netcdf_model_mask":"masks_*GRIDNAME*.nc",
	"netcdf_model_mask_lake":"flake_*GRIDNAME*.nc",
	"txt_model_grid_lat":"*GRIDNAME*_lat.txt",
	"txt_model_grid_lon":"*GRIDNAME*_lon.txt",
    "tif_copernicus_90m":"Copernicus_DSM_30_*HEMISPHERE-NS**LAT2DIGITS*_00_*HEMISPHERE-EW**LON3DIGITS*_00_DEM.tif",
    "netcdf_copernicus_lowres":"Copernicus_DSM_30_*HEMISPHERE-NS**LAT2DIGITS*_00_*HEMISPHERE-EW**LON3DIGITS*_00_DEM.tif.nc",
    "netcdf_copernicus_global_0_360":"CopernicusGlobal_float32_0_360.nc"
}

files_out = {
	"img_1km_global_orog":"1km_global_orog.png",
	"img_1km_global_orog_smooth":"1km_global_orog_smooth.png",
    "img_model_grid_orog":"model_grid_orog_*GRIDNAME*.png",
    "img_model_grid_operational_orog":"model_grid_operational_orog_*GRIDNAME*.png",
    "netcdf_model_grid_operational_orog":"model_grid_operational_orog_*GRIDNAME*.nc",
    "netcdf_model_grid_ogwd_params":"model_grid_ogwd_params_*GRIDNAME*.nc",
    "netcdf_model_grid_tofd_params":"model_grid_tofd_params_*GRIDNAME*.nc",
    "netcdf_1km_grid_orog_out":"km1_global_orog.nc"
}

files_work = {
	"netcdf_1km_smooth_orog":"1km_global_orog_smooth_to_*GRIDNAME*.nc",
	"netcdf_model_grid_orog":"model_grid_orog_*GRIDNAME*.nc"
}

model_grid = {
	"GRIDNAME":"globo_KM078",
    "gridspacing_max_km":78
}

filtering_1km_ecmwf = {
	"d":1, # [km]
    "Cfilt_max_zonal":3,
    "CN_minmax":3
}

filtering_Ndx_ecmwf = {
	"d":1, # [pixel]
    "Ndx":2 # [pixel]
}


"""
USER MANUAL

### filtering_1km_ecmwf

    THIS SECTIONS CONTAINS THE PARAMETERS USED TO SMOOTH THE RAW 1KM OROGRAPHY
    TO A 1KM RESOLUTION OROGRAPHY BUT WITH OROGRAPHIC SCALES BELOW MODEL TARGET
    GRID BOX SPACING FILTERED OUT. THE FILTERING IS PERFORMED WITH THE ECMWF
    FILTER.
    
    d: the lenght of the filter edge. [km]    

    Cfilt_max_zonal: defines the maximum filterng scale in zonal direction; 
                     the number of longitudinal points / Cfilt_max_zonal is 
                     the maximum filterng scale in zonal direction. This is
                     needed in order to avoid to use too many points near the
                     poles. [dimensionless]
                     
    CN_minmax: defines the dimension of the neighbrhood considered by the 
               local_minmax algorithm; 
               the filter lenght scale D / CN_minmax is the dimension of the
               local neighbrhood. The greater CN_minmax, the greater is the 
               effect of the local_minmax algorithm and less is the effect
               of the whole filtering precedure (the final smoothed orog will
               resemble the original one) [dimensionless]

### filtering_Ndx_ecmwf

    THIS SECTIONS CONTAINS THE PARAMETERS USED TO ADDITIONALLY SMOOTH THE MODEL
    MEAN OROGRAPHY. THIS CAN BE NECESSARY FOR EXAMPLE BECAUSE OF MODEL INSTABILITIES
    WITH TOO STEEP OROGRAPHY, OR TO MATCH THE EFFECTIVE OROGRAPHIC RESOLUTION
    OF THE MODEL (KANEHAMA ET AL)

    d: the lenght of the filter edge. [pixels]    
    
    Ndx: filtering scale of the ecmwf filter; coincides with D. [pixels]
    
"""
