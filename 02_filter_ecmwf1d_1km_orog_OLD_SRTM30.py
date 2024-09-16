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
import oroglobo_plotting as oroplot
from scipy import signal
import utils.filtering_ecmwf as orofilt_ecmwf
import utils.filtering_common as orofilt_common
from haversine import haversine
import scipy.ndimage as ndimage


######################## IMPORT PARAMETERS ##########################
gridname=oropar.model_grid["GRIDNAME"]
gridspacing_max_km=oropar.model_grid["gridspacing_max_km"]


path_data_in=oropar.paths_in["srtm30_data"]
path_img_out=oropar.paths_out["img_out"].replace("*GRIDNAME*", gridname) 
path_data_out=oropar.paths_work["workdir_grid"].replace("*GRIDNAME*", gridname) 

netcdf_orog_in=oropar.files_in["netcdf_srtm30_global_nan_to_zero_0_360"]

img_1km_global_raw_out=oropar.files_out["img_1km_global_orog"]
img_1km_global_smooth_out=oropar.files_out["img_1km_global_orog_smooth"]
netcdf_1km_smooth_orog_out=oropar.files_work["netcdf_1km_smooth_orog"].replace("*GRIDNAME*", gridname) 

d=int(oropar.filtering_1km_ecmwf["d"])
Cfilt_max_zonal=int(oropar.filtering_1km_ecmwf["Cfilt_max_zonal"])
CN_minmax=int(oropar.filtering_1km_ecmwf["CN_minmax"])

#####################################################################

# Open raw 1km orography file
data_orog = xr.open_dataset(path_data_in+netcdf_orog_in)
lat = data_orog.latitude.values
lon = data_orog.longitude.values
nlat=len(lat)
nlon=len(lon)

#### DATA_OROG.ELEV: FIRST DIMENSION (AXIS 0): LAT; SECOND DIMENSION (AXIS 1): LON. 


# plot
oroplot.orography_plot(data_orog.elev,path_img_out+img_1km_global_raw_out,2400)
print("oroginal plot: done!")

#####################################################################
########################## FILTERING ################################

filt_scale_km=float(gridspacing_max_km) # maximum model target grid spacing, km

# d is read at the beginning of the script as an external parameter

#### determine d and D, converting from km to pixels
#### they are uniform along longitudes, but varies with latitude due to meridian convergence

# for longitudes:
#   1) determine the lenght in km of each pixel in meridional direction
#   2) determine how many pixel are needed to "cover" the filtering scale 
#   3) set D parameter

dlat=180/len(lat)
km_in_each_pixel_meridional=haversine((0,0),(dlat,0),unit='km') # this is the distance, in km, of two consecutive latitudes (it is uniform since the grid is regular)
filtering_scale_pixel_meridional=int(filt_scale_km/km_in_each_pixel_meridional) # the filtering scale in the meridional direction, in pixels
D_meridional_pixel=np.full(len(lon),filtering_scale_pixel_meridional) # filtering parameter D for meridional direction; is always the same

# at each latitude:
#   1) determine the lenght in km of each pixel in zonal direction
#   2) determine how many pixel are needed to "cover" the filtering scale 
#   3) set rc parameter

D_zonal_pixel=np.zeros(len(lat))
dlon=360/len(lon)

for i in range(len(lat)):
    
    km_in_each_pixel_zonal=haversine((lat[i],lon[0]),(lat[i],lon[1]),unit='km') # this is the distance, in km, of two consecutive longitudes, at each latitude
    filtering_scale_pixel_zonal=int(filt_scale_km/km_in_each_pixel_zonal) # the filtering scale in the zonal direction, in pixels
    
    #### DEFINE A MAXIMUM FILTERNG SCALE IN ZONAL DIRECTION (OTHWERWISE AT POLES I WOULD HAVE UNFEASIBLE SCALES, TOO BIG) 
    #### Cfilt_max_zonal is an external TUNING PARAMETER
    if filtering_scale_pixel_zonal>int(len(lon)/Cfilt_max_zonal):
        D_zonal_pixel[i]=int(len(lon)/Cfilt_max_zonal)
    else:
        D_zonal_pixel[i]=filtering_scale_pixel_zonal

print('D calculation: done!')


#### FILTERING ALONG axis 0 (at each lon, filter along MERIDIONAL DIRECTION)

orogsmooth0=np.zeros((nlat,nlon))

for j in range(nlon):
    
    D=D_meridional_pixel[j]
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    orogsmooth0[:,j]=ndimage.convolve1d(data_orog.elev[:,j],F1,mode='mirror')  # the edges are the poles --> approximate with "mirror"

print('meridional filtering: done!')

#### FILTER THE PREVIOUS FIELD ALONG AXIS 1 (at each LAT, filter along ZONAL DIRECTION)

orogsmooth=np.zeros((nlat,nlon))

for i in range(nlat):
    
    D=D_zonal_pixel[i]
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    orogsmooth[i,:]=ndimage.convolve1d(orogsmooth0[i,:],F1,mode='wrap')  # the edges are datelines or greenwich --> wrap around the globe
    
print('zonal filtering: done!')

del orogsmooth0 # not needed anymore, free some memory


#### APPLY LOCAL MINMAX CONSTRAINT
# define neighborhood lenghtscale for each lat/lon
# DEFINED AS A FRACTION OF THE FILTERING SCALE D

# CN_minmax is an external TUNING PARAMETER

N_minmax_zonal_pixel=(D_zonal_pixel/CN_minmax).astype(int)
N_minmax_meridional_pixel=(D_meridional_pixel/CN_minmax).astype(int)

orogsmooth=orofilt_common.apply_local_minmax_constraint_2D_variableN(orogsmooth, data_orog.elev, N_minmax_zonal_pixel, N_minmax_meridional_pixel, ['wrap','symmetric'])

print('local minmax constraint: done!')

######################################################################
############## PLOT SMOOTH 1KM OROG AND SAVE TO NETCDF ###############

oroplot.orography_plot(orogsmooth,path_img_out+img_1km_global_smooth_out,2400)
print('smooth plot: done!')

orogsmooth_da=xr.DataArray(orogsmooth, coords=[('latitude', lat),('longitude', lon)])
# save as netcdf
orogsmooth_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_smooth_orog_out)

print('save to netcdf: done!')


# run in approx. 15 min on tintin