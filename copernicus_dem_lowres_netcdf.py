#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:50:41 2024

@author: guidodavoli


AT PRESENT SOME OF THE ORIGINAL TILES ARE MISSING



this script takes the original copernicus 90m dem in tif format
and transofrm it in a 900m resolution dataset in netcdf format.
file naming conventions are the same, BUT:

- files filled with zeros are added if original data are not present
- all files have the same number of pixels: upsampling is performed
    for original data at lat >+50 (the scales truly represented concide
    with the original dataset, no added information)
    
- THE NETCDF FILE STRUCTURE CAN BE IMPROVED (PANOPLY DOES NOT READ THEM WELL, CDO DOES)


#### COPERNICUS90m data structure

divided in four regions: NE, NW, SE, SW

tiles names are the coordinates of the "lower-left" corner

1 tile is 1deg x 1deg

NE: tiles start from N00 .... E000; tiles ends at N89 ... E179
NW: tiles start from N00 .... W001; tiles ends at N89 ... W180
SE: tiles start from S01 .... E000; tiles ends at S90 ... E179
SW: tiles start from S01 .... W001; tiles ends at S90 ... W180


"""

import xarray as xr
import rioxarray as rioxr
import matplotlib.pyplot as plt
import glob
import numpy as np
import oroglobo_parameters as oropar


# IMPORT PARAMETERS

path_to_data_original_template=oropar.paths_in["copernicus_90m"]
path_to_data_lowres_template=oropar.paths_in["copernicus_lowres"]

files_in_original_template=oropar.files_in["tif_copernicus_90m"]


def preproc_ds(ds):
    
    ds=ds.isel( x=slice(0,len(ds.x)-1) , y=slice(0,len(ds.y)-1) ) # exclude the last lat and lon, which are repetitions of the first point in adjacent tiles.
    ds=ds.rename({'x':'longitude','y':'latitude'}) # change names
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
    """    
    else:
        
        # CAMBIA COMUNQUE LE LON CON QUELLE CALCOLATE SOPRA
        print("longitudes correction")
        
        dlat=abs(ds.latitude.values[1]-ds.latitude.values[0]) # this is the spacing in lat in deegress. this is always the same in the dataset; we want this spacing also in longitudes
        
        lonmin=ds.longitude.values.min()
        lonmax=ds.longitude.values.max()
        
        #print(lonmin,lonmax)
        newlon=np.linspace(lonmin, lonmin+1.0, nlat+1)[:-1] # trick to get the correct list of points
        
        ds['longitude']=newlon
    """
    #print(ds.longitude.values,len(ds.longitude.values))
    
    #print(lonmin)
    #print(np.max(np.abs(np.diff(ds.longitude.data))))
    
    return ds


#sectors=['NE','NW','SE','SW']
sectors=['SE','SW']

for sec in sectors:
    
    """
    NE: tiles start from N00 .... E000; tiles ends at N89 ... E179
    NW: tiles start from N00 .... W001; tiles ends at N89 ... W180
    SE: tiles start from S01 .... E000; tiles ends at S90 ... E179
    SW: tiles start from S01 .... W001; tiles ends at S90 ... W180
    
    these are the lower-left corners of the tiles
    but the pixel inside have lats and lon starting
    from the upper left corner! (lat+1)
    
    """
    
    if sec=='NE':
        lats=np.arange(0,90,1) # 0 to 89
        lons=np.arange(0,180,1) # 0 to 179
    if sec=='NW':
        lats=np.arange(0,90,1) # 0 to 89
        lons=np.arange(-1,-181,-1) # -1 to -180
    if sec=='SE':
        lats=np.arange(-1,-91,-1) # -1 to -90
        lons=np.arange(0,180,1) # 0 to 180
    if sec=='SW':
        lats=np.arange(-1,-91,-1) # -1 to -90
        lons=np.arange(-1,-181,-1) # -1 to -180

    path_to_data_original=path_to_data_original_template.replace("*SECTOR*", sec) 
    path_to_data_lowres=path_to_data_lowres_template.replace("*SECTOR*", sec) 

    files_list_original = glob.glob(path_to_data_original+'Copernicus_DSM_*.tif') # all the files for this hemisphere

    for lat in lats:
        for lon in lons:
            
            tif_tile_name=files_in_original_template
            tif_tile_name=tif_tile_name.replace("*HEMISPHERE-NS*", sec[0])
            tif_tile_name=tif_tile_name.replace("*HEMISPHERE-EW*", sec[1])
            tif_tile_name=tif_tile_name.replace("*LAT2DIGITS*", str(abs(lat)).zfill(2))
            tif_tile_name=tif_tile_name.replace("*LON3DIGITS*", str(abs(lon)).zfill(3))
            
            file_original=path_to_data_original+tif_tile_name
            
            #print(file_original)
            
            if file_original in files_list_original: # the tile exists --> land (at least partially)
                
                filetif=file_original
                #ds=xr.open_mfdataset(files_list, combine="by_coords", engine='rasterio', preprocess=preproc_ds, parallel=True) # is parallel=True necessary?
                tif_data = preproc_ds( rioxr.open_rasterio(filetif).to_dataset(name='elev'))
                
                #tif_data.to_zarr(filetif+'.zarr',mode='a')
                
                # UPSCALING DATA TO 900M RESOLUTION (10 TIMES)
                llats=tif_data.latitude.data[::12].astype('float32')
                llons=tif_data.longitude.data[::12].astype('float32')
                tif_data_lowres=tif_data.interp(latitude=llats,longitude=llons,method='linear')
                
                #tif_data_lowres.to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float16', "zlib": True, "complevel": 5}})
                #tif_data_lowres.to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {"zlib": True, "complevel": 5}})
                
                # instead of writing directly to netcdf, pas through dataarray. Original files labels, info, bands... are lost, but necdf file are lighter and more readable.
                tif_data_lowres_da=xr.DataArray(tif_data_lowres.elev.data[0], coords=[('latitude', tif_data_lowres.latitude.data),('longitude', tif_data_lowres.longitude.data)])
                tif_data_lowres_da.to_dataset(name = 'elev').to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float16', "zlib": True, "complevel": 5}})
                #tif_data_lowres_da.to_dataset(name = 'elev').to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {"zlib": True, "complevel": 5}})
                
                # engine h5netcdf permette di usare float16, zlib, compression; è più lento ma genera file + piccoli
                # poi cdo collgrid mi genera comunque un file grosso come se non avessi usato nessuna compressione
                # senza compressione è:
                #     tif_data.to_netcdf('netcdf/'+filetif.replace(path_to_data, '')+'.nc',encoding={'elev': {'dtype': 'float32'}})
                ### FUNZIONA, POI DA FUORI CDO COLLGRID --> FUNZIONA SE NON HO BUCHI NEI TILES (ZONE DI MARE DEVO CREARE NC FITTIZI CON ZERI)
                
                # con cdo -z zip,9 compimo anche file finale!!! cdo -z zip,9 -collgrid Copernicus_DSM_30_N* megafile.nc
                # con zip,9 panoply inifnitamente lento ad aprire
                # con zip,5 dimensone quasi uguale e tempi molto più rapidi
                
                # senza dtype: float16 non mi da problemi con panoply; con, panoply non legge
                # BUT DTYPE FLOAT16 IS VERY IMPORTANT TO TRUNCATE LAT AND LONS AND TO AVOID TRUCATION ERRORS IN COORDINATES WHEN USING THE GENERATED FILES LATER

            else: # the tile not exists --> ocean --> generate a .nc filled with zeros            
            
                filetif=file_original
                tif_data_example = preproc_ds( rioxr.open_rasterio(files_list_original[0]).to_dataset(name='elev')) # open existing tif file as example
                
                nlat_netcdf=len(tif_data_example.latitude.data[::12])
                nlon_netcdf=len(tif_data_example.longitude.data[::12])
                
                if nlon_netcdf!=nlat_netcdf:
                    
                    print("SOMETHING WRONG")
                    
                lat_netcdf=np.linspace(lat+1.0, lat, nlat_netcdf+1)[:-1] # trick to get the correct list of points
                lon_netcdf=np.linspace(lon, lon+1.0, nlon_netcdf+1)[:-1] # trick to get the correct list of points
                         
                lat_netcdf=lat_netcdf.astype('float32')
                lon_netcdf=lon_netcdf.astype('float32')
                
                datazeros=np.zeros((nlat_netcdf,nlon_netcdf))
                
                tif_data_lowres_da=xr.DataArray(datazeros, coords=[('latitude', lat_netcdf),('longitude', lon_netcdf)])
                tif_data_lowres_da.to_dataset(name = 'elev').to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float16', "zlib": True, "complevel": 5}})
                #tif_data_lowres_da.to_dataset(name = 'elev').to_netcdf(path_to_data_lowres+filetif.replace(path_to_data_original, '')+'.nc',engine="h5netcdf", encoding={'elev': {"zlib": True, "complevel": 5}})
                
                
            #print(np.max(np.abs(np.diff(tif_data_lowres_da.longitude.data))))
            