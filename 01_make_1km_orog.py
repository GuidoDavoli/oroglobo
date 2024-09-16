#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:12:51 2023

@author: guidodavoli

THIS PROGRAM:

1) LOAD THE APPROX. 1km "LOW RESOLUTION" COPERNICUS DATA WITH XARRAY
2) CONVERTS LONGITUDES TO 0..360
3) SAVE TO A UNIQUE NETCDF THE GLOABL 1KM OROGRAPHY


#### COPERNICUS90m data structure

divided in four regions: NE, NW, SE, SW

tiles names are the coordinates of the "lower-left" corner

1 tile is 1deg x 1deg

NE: tiles start from N00 .... E000; tiles ends at N89 ... E179
NW: tiles start from N00 .... W001; tiles ends at N89 ... W180
SE: tiles start from S01 .... E000; tiles ends at S90 ... E179
SW: tiles start from S01 .... W001; tiles ends at S90 ... W180

EXAMPLES: 
    
x=LatLon23.string2latlon('45 S','25 W','d% %H')
x.lon.range360()

y=LatLon23.LatLon(lat=-30,lon=-20)
y.to_string('d% %H')

p1 = Polygon([(0,0), (1,1), (1,0)])
p2 = Polygon([(0,1), (1,0), (1,1)])
print(p1.intersects(p2))

THE DEM TILES CORRESPONDS TO A 1X1 DEG GRID
--> CHECK INTERSECTIONS WITH THESE SQUARES (0..360 FORMAT?) 
AND THEN CONVERT TO HEMISPHERE FORMAT AND OPEN FILE

TREAT THE GLOBE AS A 2D GRID, ANCHE A SX DI -180 E A DX DI +180

"""

import numpy as np
import xarray as xr
import oroglobo_parameters as oropar
import oroglobo_functions as orofunc
import oroglobo_plotting as oroplot
from shapely.geometry import Polygon
import glob
import time


def preproc_ds(ds):
    
    #print(ds)
    #ds=ds.isel( x=slice(0,len(ds.x)-1) , y=slice(0,len(ds.y)-1) ) # exclude the last lat and lon, which are repetitions of the first point in adjacent tiles.
    #ds=ds.rename({'x':'longitude','y':'latitude'}) # change names
    #ds=ds.rename({'band_data':'elev'}) # change names
    
    nlat=len(ds.latitude)
    nlon=len(ds.longitude)
    
    if nlon<nlat:
        # it means that we are above 50N/50S --> non uniform resolution.
        print("!!!!!!!!!!!!!!!!")
        
        dlat=abs(ds.latitude.values[1]-ds.latitude.values[0]) # this is the spacing in lat in deegress. this is always the same in the dataset; we want this spacing also in longitudes
        
        lonmin=ds.longitude.values.min()
        lonmax=ds.longitude.values.max()
        
        #print(lonmin,lonmax)
        newlon=np.linspace(lonmin, lonmin+1.0, nlat+1)[:-1] # trick to get the correct list of points
    
        #print(newlon,len(newlon))
        
        ds=ds.interp(longitude=newlon,method='linear',kwargs={"fill_value":"extrapolate"})
    
    #print(ds.longitude.values,len(ds.longitude.values))
    
    return ds


def wrap360(ds, lon='lon'):
    """
    
    https://github.com/pydata/xarray/issues/577
    
    wrap longitude coordinates to 0..360

    Parameters
    ----------
    ds : Dataset
        object with longitude coordinates
    lon : string
        name of the longitude ('lon', 'longitude', ...)

    Returns
    -------
    wrapped : Dataset
        Another dataset array wrapped around.
    """

    # wrap -180..179 to 0..359    
    ds.coords[lon] = np.mod(ds[lon], 360)

    # sort the data
    return ds.reindex({ lon : np.sort(ds[lon])})



# IMPORT PARAMETERS

path_data_in=oropar.paths_in["copernicus_lowres"]
path_data_out=oropar.paths_out["1km_out"]

path_img_out=oropar.paths_out["1km_out"]

netcdf_1km_grid_orog_out=oropar.files_out["netcdf_1km_grid_orog_out"]
img_1km_global_orog_out=oropar.files_out["img_1km_global_orog"]
netcdf_copernicus_lowres_in=oropar.files_in["netcdf_copernicus_lowres"]


##### LOAD DATA WITH XARRAY

sectors=['NE','NW','SE','SW']

files_list_lowres=[]

for sec in sectors:

    path_to_data_lowres=path_data_in.replace("*SECTOR*", sec) 

    files_list_lowres.extend(glob.glob(path_to_data_lowres+'Copernicus_DSM_*.tif.nc')) # all the files for this sector


start = time.time()

ds=xr.open_mfdataset(files_list_lowres, combine="by_coords", preprocess=preproc_ds, parallel=True)

end = time.time()
print("Open_mfdataset: ",str(end - start), ' s')

#dataplot=ds.sel( longitude=slice(10.5,11.7) , latitude=slice(49.2,47.4) ).fillna(0).elev.data[0]
#oroplot.orography_plot(dataplot, 'figprova.png', 400)

#### convert longitudes to 0...360

ds=wrap360(ds,'longitude')

#start = time.time()

# così pesa meno ma panoply non legge
#ds.to_netcdf(path_data_out+'CopernicusGlobal.nc',engine="h5netcdf", encoding={'elev': {"zlib": True, "complevel": 5}})
# così pesa di più ma panoply legge
ds.to_netcdf(path_data_out+'CopernicusGlobal_float32_0_360.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float32', "zlib": True, "complevel": 5}})

#end = time.time()
#print("save netcdf: ",str(end - start), ' s')

# TRY: DIFFERENT CHUNKSIZE FOR LOADING DATA TO SEE IF THERE IS A SPEEDUP;
#      DIFFERENT OPTIONS FOR NETCDF SAVING TO SEE IF THERE IS A SPEEDUP/SPACE SAVING
# tried some combinations but automatic parameters are the best





















