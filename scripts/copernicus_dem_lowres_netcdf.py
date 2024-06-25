#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:50:41 2024

@author: guidodavoli

this script takes the original copernicus 90m dem in tif format
and transofrm it in a 900m resolution dataset in netcdf format.
file naming conventions are the same.


"""

import xarray as xr
import rioxarray as rioxr
import matplotlib.pyplot as plt
import glob
import numpy as np

def preproc_ds(ds):
    
    ds=ds.isel( x=slice(0,len(ds.x)-1) , y=slice(0,len(ds.y)-1) ) # exclude the last lat and lon, which are repetitions of the first point in adjacent tiles.
    ds=ds.rename({'x':'longitude','y':'latitude'}) # change names
    #ds=ds.rename({'band_data':'elev'}) # change names
    
# =============================================================================
#     nlat=len(ds.latitude)
#     nlon=len(ds.longitude)
#     
#     if nlon<nlat:
#         # it means that we are above 50N/50S --> non uniform resolution.
#         print("!!!!!!!!!!!!!!!!")
#         
#         dlat=abs(ds.latitude.values[1]-ds.latitude.values[0]) # this is the spacing in lat in deegress. this is always the same in the dataset; we want this spacing also in longitudes
#         
#         lonmin=ds.longitude.values.min()
#         lonmax=ds.longitude.values.max()
#         
#         #print(lonmin,lonmax)
#         newlon=np.linspace(lonmin, lonmin+1.0, nlat+1)[:-1] # trick to get the correct list of points
#     
#         #print(newlon,len(newlon))
#         
#         ds=ds.interp(longitude=newlon,method='linear',kwargs={"fill_value":"extrapolate"})
#     
#     #print(ds.longitude.values,len(ds.longitude.values))
# =============================================================================
    
    return ds


path_to_data_original_template='/work_big/users/davoli/copernicus_dem_90m_*HEMISPHERE*/'
path_to_data_lowres_template='/work_big/users/davoli/copernicus_dem_lowres/copernicus_dem_90m_*HEMISPHERE*/'

hemispheres=['NE','NW','SE','SW']

for hem in hemispheres:

    path_to_data_original=path_to_data_original_template.replace("*HEMISPHERE*", hem) 
    path_to_data_lowres=path_to_data_lowres_template.replace("*HEMISPHERE*", hem) 

    files_list_original = glob.glob(path_to_data_original+'Copernicus_DSM_*.tif') # all the files for this hemisphere


    for filetif in files_list_original:
        
        #ds=xr.open_mfdataset(files_list, combine="by_coords", engine='rasterio', preprocess=preproc_ds, parallel=True) # is parallel=True necessary?
        tif_data = preproc_ds( rioxr.open_rasterio(filetif).to_dataset(name='elev'))
        
        #tif_data.to_zarr(filetif+'.zarr',mode='a')
        
        # UPSCALING DATA TO 900M RESOLUTION (10 TIMES)
        tif_data_lowres=tif_data.interp(latitude=tif_data.latitude.data[::10],longitude=tif_data.longitude.data[::10],method='linear')
        
        tif_data_lowres.to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float16', "zlib": True, "complevel": 5}})
        
        # engine h5netcdf permette di usare float16, zlib, compression; è più lento ma genera file + piccoli
        # poi cdo collgrid mi genera comunque un file grosso come se non avessi usato nessuna compressione
        # senza compressione è:
        #     tif_data.to_netcdf('netcdf/'+filetif.replace(path_to_data, '')+'.nc',encoding={'elev': {'dtype': 'float32'}})
        ### FUNZIONA, POI DA FUORI CDO COLLGRID --> FUNZIONA SE NON HO BUCHI NEI TILES (ZONE DI MARE DEVO CREARE NC FITTIZI CON ZERI)
        
        # con cdo -z zip,9 compimo anche file finale!!! cdo -z zip,9 -collgrid Copernicus_DSM_30_N* megafile.nc
        # con zip,9 panoply inifnitamente lento ad aprire
        # con zip,5 dimensone quasi uguale e tempi molto più rapidi
        