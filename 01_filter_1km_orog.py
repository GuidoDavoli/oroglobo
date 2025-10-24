#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

This code:
    
    - filters the Copernicus (approx.) 1km global orography dataset
    to a scale corresponding to the maximum grid spacing S of the 
    target operational model grid. The result is a (approx.) 1km
    resolution map, but with all scales below S filtered out.
    Creating this field allows to filter uniformly over the globe 
    the orography to a resolution consistent with the maximum 
    grid spacing S of the model target grid; in this way, 
    in the next step we will obtain a mean orography which represents
    the same scales over the globe even if the lat-lon grid spacing
    is not uniform

"""

import numpy as np
import xarray as xr
import oroglobo_plotting as oroplot
import utils.filtering_ecmwf as orofilt_ecmwf
import utils.filtering_common as orofilt_common
from haversine import haversine
import scipy.ndimage as ndimage
from multiprocessing import Pool
import time
import yaml


######################## IMPORT PARAMETERS ##########################

configname='oroglobo_parameters.yaml'
with open(configname, 'r', encoding='utf-8') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

gridname=cfg['model_grid']['GRIDNAME']
gridspacing_max_km=cfg['model_grid']['gridspacing_max_km']


path_data_in=cfg['paths_in']['copernicus_lowresnc_global']
path_img_out=cfg['paths_out']['img_out'].replace("*GRIDNAME*", gridname) 
path_data_out=cfg['paths_work']['workdir_grid'].replace("*GRIDNAME*", gridname) 

netcdf_orog_in=cfg['files_in']['netcdf_copernicus_lowres_global']

img_1km_global_raw_out=cfg['files_out']['img_1km_global_orog']
img_1km_global_smooth_out=cfg['files_out']['img_1km_global_orog_smooth']
netcdf_1km_smooth_orog_out=cfg['files_work']['netcdf_1km_smooth_orog'].replace("*GRIDNAME*", gridname) 

d=int(cfg['filtering_1km_ecmwf']['d'])
Cfilt_max_zonal=int(cfg['filtering_1km_ecmwf']['Cfilt_max_zonal'])
CN_minmax=int(cfg['filtering_1km_ecmwf']['CN_minmax'])

Nparal=int(cfg['parallel_execution']['Nparal'])
img_dpi=int(cfg['plotting']['dpi'])
make_plot=bool(cfg['plotting']['make_plot'])


#####################################################################

start = time.time()

# Open raw 1km orography file
data_orog = xr.open_dataset(path_data_in+netcdf_orog_in)
lat = data_orog.latitude.values
lon = data_orog.longitude.values
nlat=len(lat)
nlon=len(lon)

#### DATA_OROG.ELEV: FIRST DIMENSION (AXIS 0): LAT; SECOND DIMENSION (AXIS 1): LON. 

# plot
if make_plot:
    oroplot.orography_plot(data_orog.elev,path_img_out+img_1km_global_raw_out,img_dpi)
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

print('D calculation...')

for i in range(len(lat)):
    
    km_in_each_pixel_zonal=haversine((lat[i],lon[0]),(lat[i],lon[1]),unit='km') # this is the distance, in km, of two consecutive longitudes, at each latitude
    filtering_scale_pixel_zonal=int(filt_scale_km/km_in_each_pixel_zonal) # the filtering scale in the zonal direction, in pixels
    
    #### DEFINE A MAXIMUM FILTERNG SCALE IN ZONAL DIRECTION (avoid problems at poles) 
    #### Cfilt_max_zonal is an external parameter
    if filtering_scale_pixel_zonal>int(len(lon)/Cfilt_max_zonal):
        D_zonal_pixel[i]=int(len(lon)/Cfilt_max_zonal)
    else:
        D_zonal_pixel[i]=filtering_scale_pixel_zonal

print('... done!')
end = time.time()
print("time: {}".format(end-start)," s\n")


#### FILTERING ALONG axis 0 (at each lon, filter along MERIDIONAL DIRECTION)


print('meridional filtering...')

orogsmooth0=np.zeros((nlat,nlon)).astype('int16')

def funz1(j):
    
    D=D_meridional_pixel[j]
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    orogsmooth0_j=ndimage.convolve1d(data_orog.elev[:,j],F1,mode='mirror')  # the edges are the poles --> approximate with "mirror"
    #print(lon[j])
    return orogsmooth0_j


pool = Pool(processes=Nparal)

####### PARALLEL COMPUTATION #######
results = [pool.apply_async(funz1, [val]) for val in range(nlon)]
for idx, val in enumerate(results):
    orogsmooth0[:,idx]= val.get()
#######

pool.close()

print('... done!')
end = time.time()
print("time: {}".format(end-start)," s\n")


#### FILTER THE PREVIOUS FIELD ALONG AXIS 1 (at each LAT, filter along ZONAL DIRECTION)
print('zonal filtering...')

orogsmooth=np.zeros((nlat,nlon)).astype('int16')

for i in range(nlat):
    
    D=D_zonal_pixel[i]
    F1=orofilt_ecmwf.filter_ECMWF_1D(d,D)
    orogsmooth[i,:]=ndimage.convolve1d(orogsmooth0[i,:],F1,mode='wrap')  # the edges are datelines or greenwich --> wrap around the globe
    
print('... done!')
end = time.time()
print("time: {}".format(end-start)," s\n")

del orogsmooth0 # not needed anymore, free some memory


#### APPLY LOCAL MINMAX CONSTRAINT
# define neighborhood lenghtscale for each lat/lon
# DEFINED AS A FRACTION OF THE FILTERING SCALE D

# CN_minmax is an external parameter

print('local minmax constrain...')

N_minmax_zonal_pixel=(D_zonal_pixel/CN_minmax).astype(int)
N_minmax_meridional_pixel=(D_meridional_pixel/CN_minmax).astype(int)

orogsmooth=orofilt_common.apply_local_minmax_constraint_2D_variableN(orogsmooth, data_orog.elev, N_minmax_zonal_pixel, N_minmax_meridional_pixel, ['wrap','symmetric'], Nparal)

print('... done!')
end = time.time()
print("time: {}".format(end-start)," s\n")

######################################################################
############## PLOT SMOOTH 1KM OROG AND SAVE TO NETCDF ###############

if make_plot:
    oroplot.orography_plot(orogsmooth,path_img_out+img_1km_global_smooth_out,img_dpi)
    print('smooth plot: done!')

orogsmooth_da=xr.DataArray(orogsmooth, coords=[('latitude', lat),('longitude', lon)])

# reverse the order of latitudes (compatibility with following scripts)
orogsmooth_da=orogsmooth_da.isel(latitude=slice(None, None, -1))

# save as netcdf
print('save to netcdf...')

orogsmooth_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_smooth_orog_out,engine="h5netcdf", encoding={'elev': {'dtype': 'int16', "zlib": True, "complevel": 5}})

print('... done!')

print('END!')

end = time.time()
print("Total time: {}".format(end-start)," s\n")
