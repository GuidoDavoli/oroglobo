#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:12:51 2023

@author: guidodavoli


1) LOAD MODEL GRID; CYCLE OVER MODEL GRIDBOXES:
2) IDENTIFY WHICH COPERNICUS DEM TILES ARE INSIDE/INTERSECTED BY 
    THE CONSIDERED MODEL GRID
3) SELECT FROM THE SELECTED TILES ONLY THE POINTS REALLY INSIDE THE MODEL GRIDBOX
    3.1) (AT PRESENT NOT DONE) INTERPOLATE THE SELECTED 90M DATA TO A GRID WITH
    THE SAME RESOLUTION SPANNING THE MODEL GRIDBOX EDGES
4) SAMPLING THE MODEL GRIDBOX-MEAN OROGRAPHY ON THE 90M GRID (THE SAME VALUE)
5) TAKE THE DIFFERENCE BETWEEN THE SELECTED 90M DATA AND THE GRIDBOX MEAN 
    OROGRAPHY TO OBTAIN THE SUBGRID SCALE OROGRAPHY
6) DEPENDING ON THE SCHEME (VAN NIEKERK, FLOW BLOCK, ECC): FILTER OUT SCALES
    BELOW 5 KM (OPTIMIZE: SELECT A FRAME OF 90M DATA AROUND THE MODEL GRIDBOX TO
    AVOID SPECTRAL ALIASING)
7) COMPUTE THE VAN NIEKERK Fn PARAMETERS AND THE FOUR LOTT&MILLER OROG PARAMS


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
np.seterr(divide='ignore', invalid='ignore')
import xarray as xr
import rioxarray as rioxr
import oroglobo_parameters as oropar
import oroglobo_functions as orofunc
from shapely.geometry import Polygon
import LatLon23



# IMPORT PARAMETERS
gridname=oropar.model_grid["GRIDNAME"]

path_grid_in=oropar.paths_in["model_grids"]
path_data_in=oropar.paths_in["copernicus_90m"]
path_data_operational_orog_on_model_grid_in=oropar.paths_out["data_out"].replace("*GRIDNAME*", gridname) 
path_data_out=oropar.paths_out["data_out"].replace("*GRIDNAME*", gridname) 
path_img_out=oropar.paths_out["img_out"].replace("*GRIDNAME*", gridname) 

txt_model_grid_lat_in=oropar.files_in["txt_model_grid_lat"].replace("*GRIDNAME*", gridname) 
txt_model_grid_lon_in=oropar.files_in["txt_model_grid_lon"].replace("*GRIDNAME*", gridname) 

netcdf_model_grid_ogwd_params_out=oropar.files_out["netcdf_model_grid_ogwd_params"].replace("*GRIDNAME*", gridname) 
netcdf_model_grid_tofd_params_out=oropar.files_out["netcdf_model_grid_tofd_params"].replace("*GRIDNAME*", gridname) 
netcdf_model_grid_operational_orog_in=oropar.files_out["netcdf_model_grid_operational_orog"].replace("*GRIDNAME*", gridname) 
tif_copernicus_90m_in=oropar.files_in["tif_copernicus_90m"]
#img_model_grid_orog_out=oropar.files_out["img_model_grid_orog"].replace("*GRIDNAME*", gridname) 


# Open operational orography on model grid and load coordinates
operational_orog_on_model_grid_da = xr.open_dataset(path_data_operational_orog_on_model_grid_in+netcdf_model_grid_operational_orog_in)
lat_model = operational_orog_on_model_grid_da.latitude.values
lon_model = operational_orog_on_model_grid_da.longitude.values
nlat_model=len(lat_model)
nlon_model=len(lon_model)


# initialize global arrays of ogwd parameters
ogwd_F1_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_F2_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_F3_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_hamp_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_stddev_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_anisotropy_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_orientation_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
ogwd_slope_on_model_grid=np.zeros((len(lat_model),len(lon_model)))


# initialize global arrays of tofd parameters
tofd_stddev_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
tofd_anisotropy_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
tofd_orientation_on_model_grid=np.zeros((len(lat_model),len(lon_model)))
tofd_slope_on_model_grid=np.zeros((len(lat_model),len(lon_model)))


#### TRANSFORM LONGITUDES FROM 0..360 TO -180--180 ########

lon_model=lon_model-180

###########################################################


for latindex in range(nlat_model):
    
    #print(f'Time: {time.time() - start}')
    
    # SOUTH POLE POINT
    if lat_model[latindex]==-90:
        
        print("South Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex+1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        northboundary=-90+dr # lat boundary of the polar cap
        
        # collect all the points from copernicus tiles falling in the polar cap
        # usign shapely
                
        southboundary=-90
        westboundary=-180
        eastboundary=180
        
        # once the boundaries of the model grid cell has been determined,
        # determine which dem tiles to open.
        # the DEM tiles are on a 1x1 deg regular grid
        model_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
        tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon)
        print(tif_filenames_list,"SONO IO 1")
    
    #print(f'Time: {time.time() - start}')
    
    # NORTH POLE POINT    
    if lat_model[latindex]==90:
        
        print("North Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex-1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        southboundary=90-dr # lat boundary of the polar cap
        
        # collect all the points from copernicus tiles falling in the polar cap
        # usign shapely
                    
        northboundary=90
        westboundary=-180
        eastboundary=180
        
        # once the boundaries of the model grid cell has been determined,
        # determine which dem tiles to open.
        # the DEM tiles are on a 1x1 deg regular grid
        model_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
        tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon)
        print(tif_filenames_list,"SONO IO 2")

    #print(f'Time: {time.time() - start}')
    
    # THESE ARE "NORMAL" LAT-LON POINTS (NOT ON THE POLE).
    if lat_model[latindex]>-90 and lat_model[latindex]<90:
        
        #print("Latitude: ",lat_model[latindex])
        
        for lonindex in range(nlon_model):
            
            print("Latitude: ",lat_model[latindex]," Longitude: ",lon_model[lonindex])
            
            dlat_north=abs(lat_model[latindex]-lat_model[latindex+1])/2
            dlat_south=abs(lat_model[latindex]-lat_model[latindex-1])/2
            
            #print(f'Time: {time.time() - start}')
            
            # FIRST LONGITUDE POINT (IN A "-180/180" GRID)
            if lonindex==0:
                
                # if i am at the first model grid longitude, the dlon_west is this lon +360
                # minus the last lon in the array. I AM ASSUMING LONGITUDES
                # FROM -180 TO 180 (bit works the same for 0-360 grids)
                dlon_west=abs(lon_model[lonindex]+360-lon_model[-1])/2 
                dlon_east=abs(lon_model[lonindex]-lon_model[lonindex+1])/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                
                eastboundary =lon_model[lonindex]+dlon_east
                
                # check if the gridbox west boundary falls in the other hemisphere or not
                if lon_model[lonindex]-dlon_west>=-180:
                    
                    # the grid cell west border is at lon >=-180
                    westboundary =lon_model[lonindex]-dlon_west
                    #print("first cell west boundary: ",westboundary)
                    
                    # once the boundaries of the model grid cell has been determined,
                    # determine which dem tiles to open.
                    # the DEM tiles are on a 1x1 deg regular grid 
                    model_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                    tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon)
                    print(tif_filenames_list,"SONO IO 3")
                           
                else:
                    
                    # the grid cell west border goes in the western hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    westboundary=180+(lon_model[lonindex]-dlon_west+180) # the parenthesis is a negative number
                    #print("first cell west boundary: ",westboundary)  
                    
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from eastboundary to -180
                
                    westboundary_1=-180
                    eastboundary_1=eastboundary
                    model_grid_box_polygon_1=Polygon([(westboundary_1,southboundary), (eastboundary_1,southboundary), (eastboundary_1,northboundary), (westboundary_1,northboundary)])
                    tif_filenames_list_1=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon_1)
                    
                    # the second poligon goes from the calculated westboundary, which is at 180-something, up to eastboundary=+180
                
                    westboundary_2=westboundary
                    eastboundary_2=180
                    model_grid_box_polygon_2=Polygon([(westboundary_2,southboundary), (eastboundary_2,southboundary), (eastboundary_2,northboundary), (westboundary_2,northboundary)])
                    tif_filenames_list_2=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon_2)
                    
                    tif_filenames_list=tif_filenames_list_1+tif_filenames_list_2
                    print(tif_filenames_list,"SONO IO 4")
                    
                    
            #print(f'Time: {time.time() - start}')
            
            # LAST LONGITUDE POINT (IN A "-180/180" GRID)
            if lonindex==(nlon_model-1):
                
                # if i am at the last model grid lon, the dlon_east is this lon minus
                # the first lon in the array +360. I AM ASSUMING LONGITUDES
                # FROM -180 TO 180 (bit works the same for 0-360 grids)
                dlon_west=abs(lon_model[lonindex]-lon_model[lonindex-1])/2 
                dlon_east=abs(lon_model[lonindex]-(lon_model[0]+360))/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                
                westboundary =lon_model[lonindex]-dlon_west
                
                # check if the gridbox east boundary falls in the other hemisphere or not
                if lon_model[lonindex]+dlon_east<=180:
                    
                    # the grid cell east border is at lon <=180
                    eastboundary =lon_model[lonindex]+dlon_east
                    #print("last cell east boundary: ",eastboundary)
                    
                    # once the boundaries of the model grid cell has been determined,
                    # determine which dem tiles to open.
                    # the DEM tiles are on a 1x1 deg regular grid 
                    model_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                    tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon)
                    print(tif_filenames_list,"SONO IO 5")
                    
                else:
                    
                    # the grid cell west border goes in the western hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    eastboundary=(lon_model[lonindex]+dlon_east)-360 
                    #print("last cell east boundary: ",eastboundary)
                    
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from westboundary to 180
                    
                    westboundary_1=westboundary
                    eastboundary_1=180
                    model_grid_box_polygon_1=Polygon([(westboundary_1,southboundary), (eastboundary_1,southboundary), (eastboundary_1,northboundary), (westboundary_1,northboundary)])
                    tif_filenames_list_1=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon_1)
                    
                    # the second poligon goes from -180 to the calculated eastboundary, which is at -180+something
                
                    westboundary_2=-180
                    eastboundary_2=eastboundary
                    
                    model_grid_box_polygon_2=Polygon([(westboundary_2,southboundary), (eastboundary_2,southboundary), (eastboundary_2,northboundary), (westboundary_2,northboundary)])
                    tif_filenames_list_2=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon_2)
                    
                    tif_filenames_list=tif_filenames_list_1+tif_filenames_list_2
                    print(tif_filenames_list,"SONO IO 6")
            
            #print(f'Time: {time.time() - start}')
            
            # "CENTRAL" (NOT FIRST, NOT LAST) LONGITUDE POINT (IN A "-180/+180" GRID) ("EASY" POINTS)
            if lonindex>0 and lonindex<nlon_model-1: 
                
                #print(f'Time: {time.time() - start} | 1')
                
                # i am not at the first or last model grid longitude
            
                dlon_west=abs(lon_model[lonindex]-lon_model[lonindex-1])/2 
                dlon_east=abs(lon_model[lonindex]-lon_model[lonindex+1])/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                westboundary =lon_model[lonindex]-dlon_west
                eastboundary =lon_model[lonindex]+dlon_east
        
                # once the boundaries of the model grid cell has been determined,
                # determine which dem tiles to open.
                # the DEM tiles are on a 1x1 deg regular grid 
                model_grid_box_polygon=Polygon([(westboundary,southboundary), (eastboundary,southboundary), (eastboundary,northboundary), (westboundary,northboundary)])
                #print(f'Time: {time.time() - start} | 1a')
                tif_filenames_list=orofunc.get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon)
                #print(f'Time: {time.time() - start} | 1b')
                # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
                #print(tif_filenames_list,"SONO IO 7")
                
                #print(f'Time: {time.time() - start} | 2')
                
                #######################################################################################################
                ### calculate the ogwd and tofd parameters from the selected data points inside the model grid cell ###
                #######################################################################################################
             
                ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
                ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
                tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope\
                =orofunc.calculate_ogwd_and_tofd_parameters_in_model_grid_box(model_grid_box_polygon,
                                                                     tif_filenames_list,
                                                                     operational_mean_orog_in_model_grid_box
                                                                     )
                #print(f'Time: {time.time() - start} | 3')
                # store the values of the parameters for the particular
                # grid box in the corresponding 2darray
                
                ogwd_F1_on_model_grid[latindex,lonindex]  =ogwd_F1
                ogwd_F2_on_model_grid[latindex,lonindex]  =ogwd_F2
                ogwd_F3_on_model_grid[latindex,lonindex]  =ogwd_F3
                ogwd_hamp_on_model_grid[latindex,lonindex]=ogwd_hamp
                
                ogwd_stddev_on_model_grid[latindex,lonindex]     =ogwd_stddev
                ogwd_anisotropy_on_model_grid[latindex,lonindex] =ogwd_anisotropy
                ogwd_orientation_on_model_grid[latindex,lonindex]=ogwd_orientation
                ogwd_slope_on_model_grid[latindex,lonindex]      =ogwd_slope
                
                tofd_stddev_on_model_grid[latindex,lonindex]     =tofd_stddev
                tofd_anisotropy_on_model_grid[latindex,lonindex] =tofd_anisotropy
                tofd_orientation_on_model_grid[latindex,lonindex]=tofd_orientation
                tofd_slope_on_model_grid[latindex,lonindex]      =tofd_slope
    
                #print(f'Time: {time.time() - start} | 4')
    
# at the end, when all 2d arrays of parameters are filled, save netcdf on disk

##### CHECK IF THERE IS THE NEED TO CHANGE -180...180 TO 0...360

### create data arrays
ogwd_F1_on_model_grid_da=xr.DataArray(ogwd_F1_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_F2_on_model_grid_da=xr.DataArray(ogwd_F2_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_F3_on_model_grid_da=xr.DataArray(ogwd_F3_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_hamp_on_model_grid_da=xr.DataArray(ogwd_hamp_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])

ogwd_stddev_on_model_grid_da=xr.DataArray(ogwd_stddev_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_anisotropy_on_model_grid_da=xr.DataArray(ogwd_anisotropy_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_orientation_on_model_grid_da=xr.DataArray(ogwd_orientation_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
ogwd_slope_on_model_grid_da=xr.DataArray(ogwd_slope_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])

tofd_stddev_on_model_grid_da=xr.DataArray(tofd_stddev_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
tofd_anisotropy_on_model_grid_da=xr.DataArray(tofd_anisotropy_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
tofd_orientation_on_model_grid_da=xr.DataArray(tofd_orientation_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])
tofd_slope_on_model_grid_da=xr.DataArray(tofd_slope_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])


### create datasets
ogwd_parameters_ds=ogwd_F1_on_model_grid_da.to_dataset(name = 'ogwd_F1')
ogwd_parameters_ds['ogwd_F2']=ogwd_F2_on_model_grid_da
ogwd_parameters_ds['ogwd_F3']=ogwd_F3_on_model_grid_da
ogwd_parameters_ds['ogwd_hamp']=ogwd_hamp_on_model_grid_da
ogwd_parameters_ds['ogwd_stddev']=ogwd_stddev_on_model_grid_da
ogwd_parameters_ds['ogwd_anisotropy']=ogwd_anisotropy_on_model_grid_da
ogwd_parameters_ds['ogwd_orientation']=ogwd_orientation_on_model_grid_da
ogwd_parameters_ds['ogwd_slope']=ogwd_slope_on_model_grid_da

tofd_parameters_ds=tofd_stddev_on_model_grid_da.to_dataset(name = 'tofd_stddev')
tofd_parameters_ds['tofd_anisotropy']=tofd_anisotropy_on_model_grid_da
tofd_parameters_ds['tofd_orientation']=tofd_orientation_on_model_grid_da
tofd_parameters_ds['tofd_slope']=tofd_slope_on_model_grid_da

### save to netcdf
ogwd_parameters_ds.to_netcdf(path_data_out+netcdf_model_grid_ogwd_params_out)
tofd_parameters_ds.to_netcdf(path_data_out+netcdf_model_grid_tofd_params_out)