#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

This code:
    
    - load model grid and cycle over model gridboxes;
    - for each gridbox:
        * loads the corresponding operational model mean orography and 
        the 1km Copernicus orography 
        * takes the difference between the two to obtain unresoved orography 
        * performs the needed filtering (selecting above/below 5 km)
        * calculates the orographic parameters needed by parameterizations
    - saves the results in netcdf files.
    
"""

import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import xarray as xr
import oroglobo_functions as orofunc
import LatLon23
import yaml


def calc_sgo_parameters(latindex):
    
    # latindex: index of the latitude dimension
        
    # initialize single lat arrays of ogwd parameters
    ogwd_F1_on_model_grid=np.zeros(len(lon_model))
    ogwd_F2_on_model_grid=np.zeros(len(lon_model))
    ogwd_F3_on_model_grid=np.zeros(len(lon_model))
    ogwd_hamp_on_model_grid=np.zeros(len(lon_model))
    ogwd_stddev_on_model_grid=np.zeros(len(lon_model))
    ogwd_anisotropy_on_model_grid=np.zeros(len(lon_model))
    ogwd_orientation_on_model_grid=np.zeros(len(lon_model))
    ogwd_slope_on_model_grid=np.zeros(len(lon_model))
    
    # initialize single lat arrays of tofd parameters
    tofd_stddev_on_model_grid=np.zeros(len(lon_model))
    tofd_anisotropy_on_model_grid=np.zeros(len(lon_model))
    tofd_orientation_on_model_grid=np.zeros(len(lon_model))
    tofd_slope_on_model_grid=np.zeros(len(lon_model))
    
    # SOUTH POLE POINT
    if lat_model[latindex]==-90:
        
        print("South Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex+1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        northboundary=-90+dr # lat boundary of the polar cap
        
        # collect all the points from copernicus tiles falling in the polar cap
                
        southboundary=-90
        westboundary=-180
        eastboundary=180
        
        # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
        # perform a mean to get the value on the pole point (i select all longitudes at the lat of the pole)
        operational_mean_orog_in_model_grid_box=np.mean(np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],method="nearest").elev.data))
    
        all_tif_data_ds=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary,eastboundary)).rename({'longitude': 'x','latitude': 'y'})
    
        #######################################################################################################
        ### calculate the ogwd and tofd parameters from the selected data points inside the model grid cell ###
        #######################################################################################################
                        
        ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
        ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
        tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope\
        =orofunc.calculate_ogwd_and_tofd_parameters_in_model_grid_box_from1kmorog(#model_grid_box_polygon,
                                                             all_tif_data_ds,
                                                             operational_mean_orog_in_model_grid_box
                                                             )
        
        # the poles are single points in reality, so the value is the same at all longitudes
        
        ogwd_F1_on_model_grid[:]  =ogwd_F1
        ogwd_F2_on_model_grid[:]  =ogwd_F2
        ogwd_F3_on_model_grid[:]  =ogwd_F3
        ogwd_hamp_on_model_grid[:]=ogwd_hamp
        
        ogwd_stddev_on_model_grid[:]     =ogwd_stddev
        ogwd_anisotropy_on_model_grid[:] =ogwd_anisotropy
        ogwd_orientation_on_model_grid[:]=ogwd_orientation
        ogwd_slope_on_model_grid[:]      =ogwd_slope
        
        tofd_stddev_on_model_grid[:]     =tofd_stddev
        tofd_anisotropy_on_model_grid[:] =tofd_anisotropy
        tofd_orientation_on_model_grid[:]=tofd_orientation
        tofd_slope_on_model_grid[:]      =tofd_slope
        
        ###########################
        ### end of calculations ###
        ###########################
    
    
    # NORTH POLE POINT    
    if lat_model[latindex]==90:
        
        print("North Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex-1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        southboundary=90-dr # lat boundary of the polar cap
        
        # collect all the points from copernicus tiles falling in the polar cap
                    
        northboundary=90
        westboundary=-180
        eastboundary=180
        
        # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
        # perform a mean to get the value on the pole point (i select all longitudes at the lat of the pole)
        operational_mean_orog_in_model_grid_box=np.mean(np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],method="nearest").elev.data))

        all_tif_data_ds=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary,eastboundary)).rename({'longitude': 'x','latitude': 'y'})
        
        #######################################################################################################
        ### calculate the ogwd and tofd parameters from the selected data points inside the model grid cell ###
        #######################################################################################################
                        
        ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
        ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
        tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope\
        =orofunc.calculate_ogwd_and_tofd_parameters_in_model_grid_box_from1kmorog(#model_grid_box_polygon,
                                                             all_tif_data_ds,
                                                             operational_mean_orog_in_model_grid_box
                                                             )
        
        # the poles are single points in reality, so the value is the same at all longitudes
        
        ogwd_F1_on_model_grid[:]  =ogwd_F1
        ogwd_F2_on_model_grid[:]  =ogwd_F2
        ogwd_F3_on_model_grid[:]  =ogwd_F3
        ogwd_hamp_on_model_grid[:]=ogwd_hamp
        
        ogwd_stddev_on_model_grid[:]     =ogwd_stddev
        ogwd_anisotropy_on_model_grid[:] =ogwd_anisotropy
        ogwd_orientation_on_model_grid[:]=ogwd_orientation
        ogwd_slope_on_model_grid[:]      =ogwd_slope
        
        tofd_stddev_on_model_grid[:]     =tofd_stddev
        tofd_anisotropy_on_model_grid[:] =tofd_anisotropy
        tofd_orientation_on_model_grid[:]=tofd_orientation
        tofd_slope_on_model_grid[:]      =tofd_slope
        
        ###########################
        ### end of calculations ###
        ###########################
        
    
    # THESE ARE "NORMAL" LAT-LON POINTS (NOT ON THE POLE).
    if lat_model[latindex]>-90 and lat_model[latindex]<90:
        
        print("Working on latitude... ",lat_model[latindex])
        
        for lonindex in range(nlon_model):
            
            dlat_north=abs(lat_model[latindex]-lat_model[latindex+1])/2
            dlat_south=abs(lat_model[latindex]-lat_model[latindex-1])/2
                        
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
                    
                    # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                    # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                    operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
                     
                    all_tif_data_ds=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary,eastboundary)).rename({'longitude': 'x','latitude': 'y'})
                    
                else:
                    
                    # the grid cell west border goes in the eastern hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    westboundary=180+(lon_model[lonindex]-dlon_west+180) # the parenthesis is a negative number                    
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from eastboundary to -180
                    westboundary_1=-180
                    eastboundary_1=eastboundary

                    # the second poligon goes from the calculated westboundary, which is at 180-something, up to eastboundary=+180            
                    westboundary_2=westboundary
                    eastboundary_2=180
                    
                    # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                    # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                    operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
                    
                    # now there is the need to define "model_grid_box_polygon" 
                    # for the next computations; is a merge of the previous two polygons
                    
                    all_tif_data_ds1=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary_1,eastboundary_1)).rename({'longitude': 'x','latitude': 'y'})
                    all_tif_data_ds2=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary_2,eastboundary_2)).rename({'longitude': 'x','latitude': 'y'})
                    all_tif_data_ds=xr.merge([all_tif_data_ds1,all_tif_data_ds2])
        

            
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
                    
                    # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                    # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                    operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
                    
                    all_tif_data_ds=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary,eastboundary)).rename({'longitude': 'x','latitude': 'y'})
                    
                else:
                    
                    # the grid cell east border goes in the western hemisphere     
                    # TO DEAL WITH THIS PROBLEM, WE HAVE TO SPLIT THE MODEL GRID CELL
                    # IN TWO POLYGONS AND CHECK THE INTERSECTIONS TWO TIMES
                    
                    eastboundary=(lon_model[lonindex]+dlon_east)-360                     
                    
                    #SPLIT THE PROBLEM IN TWO PLYGONS
                    
                    # the first poligon goes from westboundary to 180
                    
                    westboundary_1=westboundary
                    eastboundary_1=180

                    # the second poligon goes from -180 to the calculated eastboundary, which is at -180+something
                
                    westboundary_2=-180
                    eastboundary_2=eastboundary
                    
                    # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                    # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                    operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
            
                    all_tif_data_ds1=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary_1,eastboundary_1)).rename({'longitude': 'x','latitude': 'y'})
                    all_tif_data_ds2=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary_2,eastboundary_2)).rename({'longitude': 'x','latitude': 'y'})
                    all_tif_data_ds=xr.merge([all_tif_data_ds1,all_tif_data_ds2])
                    
                        
            # "CENTRAL" (NOT FIRST, NOT LAST) LONGITUDE POINT (IN A "-180/+180" GRID) ("EASY" POINTS)
            if lonindex>0 and lonindex<nlon_model-1: 
                                
                # i am not at the first or last model grid longitude
            
                dlon_west=abs(lon_model[lonindex]-lon_model[lonindex-1])/2 
                dlon_east=abs(lon_model[lonindex]-lon_model[lonindex+1])/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                westboundary =lon_model[lonindex]-dlon_west
                eastboundary =lon_model[lonindex]+dlon_east
                
                # here I use LatLon23 range360() function as a workaraound to get lon in 0..360 format in order to get the correct model mean orography 
                # method=nearest allow for selection of the nearest point if the coordinates are not exact (can happen due to truncation errors).
                operational_mean_orog_in_model_grid_box=np.float32(operational_orog_on_model_grid_da.sel(latitude=lat_model[latindex],longitude=LatLon23.Longitude(lon_model[lonindex]).range360(),method="nearest").elev.data)
                            
                all_tif_data_ds=orog1kmraw.sel(latitude=slice(northboundary,southboundary),longitude=slice(westboundary,eastboundary)).rename({'longitude': 'x','latitude': 'y'})
                
                
            #######################################################################################################
            ### calculate the ogwd and tofd parameters from the selected data points inside the model grid cell ###
            #######################################################################################################
                            
            ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
            ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
            tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope\
            =orofunc.calculate_ogwd_and_tofd_parameters_in_model_grid_box_from1kmorog(#model_grid_box_polygon,
                                                                 all_tif_data_ds,
                                                                 operational_mean_orog_in_model_grid_box
                                                                 )
            
            ogwd_F1_on_model_grid[lonindex]  =ogwd_F1
            ogwd_F2_on_model_grid[lonindex]  =ogwd_F2
            ogwd_F3_on_model_grid[lonindex]  =ogwd_F3
            ogwd_hamp_on_model_grid[lonindex]=ogwd_hamp
            
            ogwd_stddev_on_model_grid[lonindex]     =ogwd_stddev
            ogwd_anisotropy_on_model_grid[lonindex] =ogwd_anisotropy
            ogwd_orientation_on_model_grid[lonindex]=ogwd_orientation
            ogwd_slope_on_model_grid[lonindex]      =ogwd_slope
            
            tofd_stddev_on_model_grid[lonindex]     =tofd_stddev
            tofd_anisotropy_on_model_grid[lonindex] =tofd_anisotropy
            tofd_orientation_on_model_grid[lonindex]=tofd_orientation
            tofd_slope_on_model_grid[lonindex]      =tofd_slope
            
            ###########################
            ### end of calculations ###
            ###########################
        
    return ogwd_F1_on_model_grid,ogwd_F2_on_model_grid,ogwd_F3_on_model_grid,ogwd_hamp_on_model_grid, \
           ogwd_stddev_on_model_grid,ogwd_anisotropy_on_model_grid,ogwd_orientation_on_model_grid,ogwd_slope_on_model_grid, \
           tofd_stddev_on_model_grid,tofd_anisotropy_on_model_grid,tofd_orientation_on_model_grid,tofd_slope_on_model_grid



############## IMPORT PARAMETERS

configname='oroglobo_parameters.yaml'
with open(configname, 'r', encoding='utf-8') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

gridname=cfg['model_grid']['GRIDNAME']

path_grid_in=cfg['paths_in']["model_grids"]
path_data_in=cfg['paths_in']["copernicus_lowresnc_global"]
path_data_operational_orog_on_model_grid_in=cfg['paths_out']["data_out"].replace("*GRIDNAME*", gridname) 
path_data_out=cfg['paths_out']["data_out"].replace("*GRIDNAME*", gridname) 
path_img_out=cfg['paths_out']["img_out"].replace("*GRIDNAME*", gridname) 

txt_model_grid_lat_in=cfg['files_in']["txt_model_grid_lat"].replace("*GRIDNAME*", gridname) 
txt_model_grid_lon_in=cfg['files_in']["txt_model_grid_lon"].replace("*GRIDNAME*", gridname) 

netcdf_model_grid_ogwd_params_out=cfg['files_out']["netcdf_model_grid_ogwd_params"].replace("*GRIDNAME*", gridname) 
netcdf_model_grid_tofd_params_out=cfg['files_out']["netcdf_model_grid_tofd_params"].replace("*GRIDNAME*", gridname) 
netcdf_model_grid_operational_orog_in=cfg['files_out']["netcdf_model_grid_operational_orog"].replace("*GRIDNAME*", gridname) 
netcdf_copernicus_lowres_global=cfg['files_in']["netcdf_copernicus_lowres_global"]

Nparal=int(cfg['parallel_execution']['Nparal'])



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


# Open raw 1km orography file
orog1kmraw = xr.open_dataset(path_data_in+netcdf_copernicus_lowres_global)

orog1kmraw.coords['longitude'] = (orog1kmraw.coords['longitude'] + 180) % 360 - 180
orog1kmraw = orog1kmraw.sortby(orog1kmraw.longitude)
orog1kmraw.load()



#### TRANSFORM LONGITUDES FROM 0..360 TO -180--180 ########
lon_model=lon_model-180

###########################################################

import time

start = time.time()
    
from multiprocessing import Pool

pool = Pool(processes=Nparal)

########## PARALLEL EXECUTION

results = [pool.apply_async(calc_sgo_parameters, [val]) for val in range(nlat_model)]
for idx, val in enumerate(results):
    ogwd_F1_on_model_grid[idx,:], \
    ogwd_F2_on_model_grid[idx,:], \
    ogwd_F3_on_model_grid[idx,:], \
    ogwd_hamp_on_model_grid[idx,:], \
    ogwd_stddev_on_model_grid[idx,:], \
    ogwd_anisotropy_on_model_grid[idx,:], \
    ogwd_orientation_on_model_grid[idx,:], \
    ogwd_slope_on_model_grid[idx,:], \
    tofd_stddev_on_model_grid[idx,:], \
    tofd_anisotropy_on_model_grid[idx,:], \
    tofd_orientation_on_model_grid[idx,:], \
    tofd_slope_on_model_grid[idx,:] \
    = val.get()
    
#######

pool.close()
    
    
# at the end, when all 2d arrays of parameters are filled, save netcdf on disk

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

end = time.time()
    
print("time: {}\n".format(end-start))