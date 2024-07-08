#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:12:51 2023

@author: guidodavoli

THIS VERSION USES XARRAY 

1) LOAD 1km GRID; CYCLE OVER GRIDBOXES:
2) IDENTIFY WHICH COPERNICUS DEM TILES ARE INSIDE/INTERSECTED BY 
    THE CONSIDERED GRID cell
3) SELECT FROM THE SELECTED TILES ONLY THE POINTS REALLY INSIDE THE GRIDBOX
4) calculate 1km mean orography


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

# =============================================================================
# 
# #########################################
# ######### Define the 1km grid ###########
# ## a regular latlon grid without poles ##
# #########################################
# 
# nlon_1km=np.int32(40000) # must be even
# nlat_1km=np.int32(20000) # must be even
# 
# dlon_1km=np.float32(360)/nlon_1km
# dlat_1km=np.float32(180)/nlat_1km
# 
# lon_pos_1km=np.zeros(np.int32(nlon_1km/2),dtype='float32')
# lat_pos_1km=np.zeros(np.int32(nlat_1km/2),dtype='float32')
# lon_neg_1km=np.zeros(np.int32(nlon_1km/2),dtype='float32')
# lat_neg_1km=np.zeros(np.int32(nlat_1km/2),dtype='float32')
# 
# dlon_half_1km=np.float32(dlon_1km/2)
# dlat_half_1km=np.float32(dlat_1km/2)
# 
# lon0_1km=dlon_half_1km
# lat0_1km=dlat_half_1km
# 
# for ilat in range(np.int32(nlat_1km/2)):
#     lat_pos_1km[ilat]=np.float32(+lat0_1km+ilat*dlat_1km)
#     lat_neg_1km[ilat]=np.float32(-lat0_1km-ilat*dlat_1km)
#     
# for ilon in range(np.int32(nlon_1km/2)):
#     lon_pos_1km[ilon]=np.float32(+lon0_1km+ilon*dlon_1km)
#     lon_neg_1km[ilon]=np.float32(-lon0_1km-ilon*dlon_1km)
# 
# lat_1km=np.append(np.flip(lat_neg_1km), lat_pos_1km) # from -90 to +90
# lon_1km=np.append(np.flip(lon_neg_1km), lon_pos_1km) # from -180 to + 180
# 
# 
# ##### INITIALIZE THE ARRAY TO STORE THE 1KM OROGRAPHY
# 
# meanorog_1km=np.zeros((nlat_1km,nlon_1km))
# 
# 
# =============================================================================
##### LOAD DATA WITH XARRAY


sectors=['NE','NW','SE','SW']
#sectors=['SW']

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

### non so perchè ci sono nan lungo longitudine -> dropna
#ds=ds.dropna('longitude')
#ds=ds.dropna('latitude')
#dataplot=ds.sel( longitude=slice(5.5,15.7) , latitude=slice(39.2,31.4) ).elev.data
#dataplot=ds.sel( longitude=slice(5.5,15.7) , latitude=slice(39.2,31.4) ).dropna('longitude').elev.data

#### convert longitudes to 0...360

ds=wrap360(ds,'longitude')


start = time.time()
# così pesa meno ma panoply non legge
#ds.to_netcdf(path_data_out+'CopernicusGlobal.nc',engine="h5netcdf", encoding={'elev': {"zlib": True, "complevel": 5}})
# così pesa di più ma panoply legge
ds.to_netcdf(path_data_out+'CopernicusGlobal_float32_0_360.nc',engine="h5netcdf", encoding={'elev': {'dtype': 'float32', "zlib": True, "complevel": 5}})

end = time.time()
print("save netcdf: ",str(end - start), ' s')

# TRY: DIFFERENT CHUNKSIZE FOR LOADING DATA TO SEE IF THERE IS A SPEEDUP;
#      DIFFERENT OPTIONS FOR NETCDF SAVING TO SEE IF THERE IS A SPEEDUP/SPACE SAVING

"""
#########################################
######## MAIN CALCULATIONS BEGIN ########
#########################################

previous_tif_filenames_list=[] # needed later to optimize access to files

for latindex in range(nlat_1km):
    
    if lat_1km[latindex]>-90 and lat_1km[latindex]<90:
    #if lat_1km[latindex]>0 and lat_1km[latindex]<10:
        
        #print("Latitude: ",lat_model[latindex])
        
        for lonindex in range(nlon_1km):
            
            print("Latitude: ",lat_1km[latindex]," Longitude: ",lon_1km[lonindex])
            
            dlat_north=dlat_half_1km
            dlat_south=dlat_half_1km
            
            #print(f'Time: {time.time() - start}')
            
            # FIRST LONGITUDE POINT (IN A "-180/180" GRID)
            if lonindex==0:
                
                dlon_west=dlon_half_1km
                dlon_east=dlon_half_1km
                
                northboundary=lat_1km[latindex]+dlat_north
                southboundary=lat_1km[latindex]-dlat_south
                
                eastboundary =lon_1km[lonindex]+dlon_east
                
                # check if the gridbox west boundary falls in the other hemisphere or not
                if lon_1km[lonindex]-dlon_west>=-180:
                    
                    # the grid cell west border is at lon >=-180
                    westboundary =lon_1km[lonindex]-dlon_west
                    #print("first cell west boundary: ",westboundary)
                    
                    # once the boundaries of the model grid cell has been determined,
                    # determine which dem tiles to open.
                    # the DEM tiles are on a 1x1 deg regular grid 
                    km1_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                    tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon)
                    print(tif_filenames_list,"SONO IO 3")
                           
                else:
                    
                    # the grid cell west border goes in the western hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    westboundary=180+(lon_1km[lonindex]-dlon_west+180) # the parenthesis is a negative number
                    #print("first cell west boundary: ",westboundary)  
                    
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from eastboundary to -180
                
                    westboundary_1=-180
                    eastboundary_1=eastboundary
                    km1_grid_box_polygon_1=Polygon([(westboundary_1,southboundary), (eastboundary_1,southboundary), (eastboundary_1,northboundary), (westboundary_1,northboundary)])
                    tif_filenames_list_1=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon_1)
                    
                    # the second poligon goes from the calculated westboundary, which is at 180-something, up to eastboundary=+180
                
                    westboundary_2=westboundary
                    eastboundary_2=180
                    km1_grid_box_polygon_2=Polygon([(westboundary_2,southboundary), (eastboundary_2,southboundary), (eastboundary_2,northboundary), (westboundary_2,northboundary)])
                    tif_filenames_list_2=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon_2)
                    
                    tif_filenames_list=tif_filenames_list_1+tif_filenames_list_2
                    print(tif_filenames_list,"SONO IO 4")
                    
                    
            #print(f'Time: {time.time() - start}')
            
            # LAST LONGITUDE POINT (IN A "-180/180" GRID)
            if lonindex==(nlon_1km-1):
                
                dlon_west=dlon_half_1km
                dlon_east=dlon_half_1km
                
                northboundary=lat_1km[latindex]+dlat_north
                southboundary=lat_1km[latindex]-dlat_south
                
                westboundary =lon_1km[lonindex]-dlon_west
                
                # check if the gridbox east boundary falls in the other hemisphere or not
                if lon_1km[lonindex]+dlon_east<=180:
                    
                    # the grid cell east border is at lon <=180
                    eastboundary =lon_1km[lonindex]+dlon_east
                    #print("last cell east boundary: ",eastboundary)
                    
                    # once the boundaries of the model grid cell has been determined,
                    # determine which dem tiles to open.
                    # the DEM tiles are on a 1x1 deg regular grid 
                    km1_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                    tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon)
                    print(tif_filenames_list,"SONO IO 5")
                    
                else:
                    
                    # the grid cell west border goes in the western hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    eastboundary=(lon_1km[lonindex]+dlon_east)-360 
                    #print("last cell east boundary: ",eastboundary)
                    
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from westboundary to 180
                    
                    westboundary_1=westboundary
                    eastboundary_1=180
                    km1_grid_box_polygon_1=Polygon([(westboundary_1,southboundary), (eastboundary_1,southboundary), (eastboundary_1,northboundary), (westboundary_1,northboundary)])
                    tif_filenames_list_1=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon_1)
                    
                    # the second poligon goes from -180 to the calculated eastboundary, which is at -180+something
                
                    westboundary_2=-180
                    eastboundary_2=eastboundary
                    
                    km1_grid_box_polygon_2=Polygon([(westboundary_2,southboundary), (eastboundary_2,southboundary), (eastboundary_2,northboundary), (westboundary_2,northboundary)])
                    tif_filenames_list_2=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon_2)
                    
                    tif_filenames_list=tif_filenames_list_1+tif_filenames_list_2
                    print(tif_filenames_list,"SONO IO 6")
            
            #print(f'Time: {time.time() - start}')
            
            # "CENTRAL" (NOT FIRST, NOT LAST) LONGITUDE POINT (IN A "-180/+180" GRID) ("EASY" POINTS)
            if lonindex>0 and lonindex<nlon_1km-1: 
                
                #print(f'Time: {time.time() - start} | 1')
                
                # i am not at the first or last 1km grid longitude
            
                dlon_west=dlon_half_1km
                dlon_east=dlon_half_1km
                
                northboundary=lat_1km[latindex]+dlat_north
                southboundary=lat_1km[latindex]-dlat_south
                westboundary =lon_1km[lonindex]-dlon_west
                eastboundary =lon_1km[lonindex]+dlon_east
        
                # once the boundaries of the model grid cell has been determined,
                # determine which dem tiles to open.
                # the DEM tiles are on a 1x1 deg regular grid 
                km1_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                #print(f'Time: {time.time() - start} | 1a')
                tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_grid_box(km1_grid_box_polygon)
                
                
                #print(f'Time: {time.time() - start} | 4')
                
                
            ###############################################################################################
            ### calculate the 1km mean orography from the selected data points inside the 1km grid cell ###
            ###############################################################################################
            
            if tif_filenames_list!=previous_tif_filenames_list:
                # tif files on disk are opened only if they are not already opened.
                # if the above condition is false, the filenames list and
                # data_exist were already opened and loaded in memory.
                print(tif_filenames_list)
                all_tif_data_ds,data_exist=orofunc.get_copernicus90m_data_xrds_from_tiles_list(tif_filenames_list)
            if data_exist:
                H=orofunc.calculate_mean_orog_in_grid_box(km1_grid_box_polygon,all_tif_data_ds)
            else:
                H=np.nan # if there are no data on disk --> the grid box is entirely on ocean
                                                     
            #print(f'Time: {time.time() - start} | 3')
            # store the values of the parameters for the particular
            # grid box in the corresponding 2darray
            
            meanorog_1km[latindex,lonindex] = H
            print("mean orog: ",str(H))
            
            ##########################
            ### end of calculation ###
            ##########################
    
            
            # to register the data already loaded and optimize access to file
            previous_tif_filenames_list=tif_filenames_list
            
            
# at the end, when all 2d arrays of parameters are filled, save netcdf on disk

##### CHECK IF THERE IS THE NEED TO CHANGE -180...180 TO 0...360



################################################################
############## PLOT  1KM OROG AND SAVE TO NETCDF ###############

oroplot.orography_plot(meanorog_1km,path_img_out+img_1km_global_orog_out,2400)
print('plot: done!')

meanorog_1km_da=xr.DataArray(meanorog_1km, coords=[('latitude', lat_1km),('longitude', lon_1km)])
# save as netcdf
meanorog_1km_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_1km_grid_orog_out)

print('save to netcdf: done!')


"""






















