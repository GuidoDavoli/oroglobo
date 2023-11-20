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
	"srtm30_data":"data_in/srtm30/"
}

paths_out = {
	"img_out":"img_out/",
    "data_out":"data_out/"
}

files_in = {
	"netcdf_in_srtm30_global_nan_to_zero":"srtm30_global_nan_to_zero.nc"
}

files_out = {
	"img_srtm30_global_nan_to_zero":"srtm30_orog.png",
	"img_srtm30_smooth":"srtm30_orog_smooth.png",
	"netcdf_srtm30_smooth":"srtm30_global_nan_to_zero_smoothed.nc"
}


