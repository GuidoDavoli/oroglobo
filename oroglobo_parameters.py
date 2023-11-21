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
	"model_masks":"data_in/model_masks/"
}

paths_out = {
	"img_out":"img_out/",
    	"data_out":"data_out/"
}

paths_work = {
	"workdir":"work/"
}

files_in = {
	"netcdf_srtm30_global_nan_to_zero_0_360":"srtm30_global_nan_to_zero_0_360.nc",
	"txt_model_grid_lat":"grids_globo_KM078_lat.txt",
	"txt_model_grid_lon":"grids_globo_KM078_lon.txt"	
}

files_out = {
	"img_srtm30_global_nan_to_zero":"srtm30_orog.png",
	"img_srtm30_smooth":"srtm30_orog_smooth.png"
}

files_work = {
	"netcdf_1km_smooth_orog":"srtm30_global_nan_to_zero_smoothed.nc",
	"netcdf_model_grid_orog":"model_grid_orog_globo_KM078.nc"
}






