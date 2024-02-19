#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:58:16 2023

@author: guido

THIS SCRIPT FILTERS A 1KM GLOBAL GRID (1KM AT EQUATOR) UNIFORMLY, TO A SCALE
CORRESPONDING TO THE MAXIMUMX GRID SPACING OF THE TARGET OPERATIONAL
MODEL GRID.
THE RESULT IS A 1KM (AT EQUATOR) RESOLUTION MAP, BUT WITH ALL SCALES BELOW
A GIVER VALUE FILTERED OUT.


"""

import numpy as np
import xarray as xr
import oroglobo_parameters as oropar
import oroglobo_functions as orofunc
import oroglobo_plotting as oroplot
import utils.filtering_zadra as orofilt_z
from scipy import signal
from haversine import haversine


# IMPORT PARAMETERS
gridname=oropar.model_grid["GRIDNAME"]
gridspacing_max_km=oropar.model_grid["gridspacing_max_km"]


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


############ FILTERING with Zadra 2018 method

filt_scale_km=float(gridspacing_max_km) # maximum model target grid spacing, km

### prepare variable rc and p. they vary with latitude to have a uniform filtering in km

# for longitudes:
#   1) determine the lenght in km of each pixel in meridional direction
#   2) determine how many pixel are needed to "cover" the filtering scale 
#   3) set rc parameter

dlat=180/len(lat)
km_in_each_pixel_meridional=haversine((0,0),(dlat,0),unit='km') # this is the distance, in km, of two consecutive latitudes (it is uniform since the grid is regular)
filtering_scale_pixel_meridional=int(filt_scale_km/km_in_each_pixel_meridional) # the filtering scale in the meridional direction, in pixels
rc_meridional_pixel=np.full(len(lon),filtering_scale_pixel_meridional) # filtering parameter rc for meridional direction; is always the same

# at each latitude:
#   1) determine the lenght in km of each pixel in zonal direction
#   2) determine how many pixel are needed to "cover" the filtering scale 
#   3) set rc parameter

rc_zonal_pixel=np.zeros(len(lat))
dlon=360/len(lon)

for i in range(len(lat)):
    
    km_in_each_pixel_zonal=haversine((lat[i],lon[0]),(lat[i],lon[1]),unit='km') # this is the distance, in km, of two consecutive longitudes, at each latitude
    filtering_scale_pixel_zonal=int(filt_scale_km/km_in_each_pixel_zonal) # the filtering scale in the zonal direction, in pixels
    
    #### DEFINE A MAXIMUM FILTERNG SCALE IN ZONAL DIRECTION (OTHWERWISE AT POLES I WOULD HAVE UNFEASIBLE SCALES, TOO BIG) 
    #### TUNING PARAMETER
    Cfilt_max_zonal=3
    if filtering_scale_pixel_zonal>int(len(lon)/Cfilt_max_zonal):
        rc_zonal_pixel[i]=int(len(lon)/Cfilt_max_zonal)
    else:
        rc_zonal_pixel[i]=filtering_scale_pixel_zonal


#### P parameter defined as a multiple of rc

p_meridional_pixel=Cfilt_max_zonal*rc_meridional_pixel
p_zonal_pixel     =Cfilt_max_zonal*rc_zonal_pixel


padp=int(max(p_meridional_pixel.max(),p_zonal_pixel.max()))  # to pad is important p, because is the larger value and defines the averagin window

rc_meridional_pixel_pad=np.pad(rc_meridional_pixel,(padp,padp)) # pad with zeros
rc_zonal_pixel_pad     =np.pad(rc_zonal_pixel,     (padp,padp)) # pad with zeros

p_meridional_pixel_pad =np.pad(p_meridional_pixel, (padp,padp)) # pad with zeros
p_zonal_pixel_pad      =np.pad(p_zonal_pixel,      (padp,padp)) # pad with zeros

### PREPARE FOR FILTERING

data_orog_elev_withframe=orofilt_z.add_frame_before_filtering(data_orog.elev, padp, ['symmetric','wrap']) # the input orog is supposed to be in lat, lon order

data_orog_elev_filtered_withframe=orofilt_z.LowPassFilter_2D_v2_UPDATE(data_orog_elev_withframe, 
                                                                       rc_meridional_pixel_pad, 
                                                                       rc_zonal_pixel_pad, 
                                                                       p_meridional_pixel_pad, 
                                                                       p_zonal_pixel_pad)

data_orog_elev_filtered=orofilt_z.remove_frame_after_filtering(data_orog_elev_filtered_withframe, padp)

### APPLY LOCAL MINMAX CONSTRAINT
### TUNING PARAMETER to define the dimension of the "local" neighborhood
L=int(min(rc_meridional_pixel.min(),rc_zonal_pixel.min())/5) 
orogsmooth=orofilt_z.apply_local_minmax_constraint(data_orog_elev_filtered,data_orog.elev,L)




"""
filt=orofunc.filter_ECMWF(1, 80)

orogsmooth = signal.convolve2d(data_orog.elev, filt, mode='same', boundary='symm').astype(np.float32) # IF BUNDARY IS SET TO WRAP, SPURIOUS OROGRAPHY AT THE POLES (np IS WRAPPED WIT sp AND VICE VERSA)
# ANOTHER SOLUTION WOULD BE TO USE "WRAP" TO GET CORRECT WRAPPING ALONG LONGITUDES, ADN CORRECT AT NP and SP: 0 at NP points, a simple average at SP points
"""

# plot
oroplot.orography_plot(orogsmooth,path_img_out+img_1km_global_smooth_out,2400)

#orogsmooth_da=xr.DataArray(orogsmooth.transpose, coords=data_orog.coords, dims=data_orog.dims, attrs=data_orog.attrs)
orogsmooth_da=xr.DataArray(orogsmooth, coords=[('latitude', lat),('longitude', lon)])

# save as netcdf
orogsmooth_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_smooth_orog_out)


# su tintin gira in 4 ore

## WARNING: AT PRESENT THE SMOOTHING OPERATION CREATE SPURIOUS OROGRAPHY OVER NORTH POLE