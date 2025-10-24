#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

This code:

    - takes the 1km reoslution orography smoothed to target model grid scale
    and aggregate to target model grid boxes (simple mean)


"""

import numpy as np
import xarray as xr
import oroglobo_plotting as oroplot
import yaml

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


############# IMPORT PARAMETERS
configname='oroglobo_parameters.yaml'
with open(configname, 'r', encoding='utf-8') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)
    
gridname=cfg['model_grid']['GRIDNAME']

path_grid_in=cfg['paths_in']["model_grids"]
path_data_in=cfg['paths_work']["workdir_grid"].replace("*GRIDNAME*", gridname) 
path_data_out=cfg['paths_work']["workdir_grid"].replace("*GRIDNAME*", gridname) 
path_img_out=cfg['paths_out']["img_out"].replace("*GRIDNAME*", gridname) 

txt_model_grid_lat_in=cfg['files_in']["txt_model_grid_lat"].replace("*GRIDNAME*", gridname) 
txt_model_grid_lon_in=cfg['files_in']["txt_model_grid_lon"].replace("*GRIDNAME*", gridname) 

netcdf_model_grid_orog_out=cfg['files_work']["netcdf_model_grid_orog"].replace("*GRIDNAME*", gridname) 
netcdf_1km_smooth_orog_in=cfg['files_work']["netcdf_1km_smooth_orog"].replace("*GRIDNAME*", gridname) 
img_model_grid_orog_out=cfg['files_out']["img_model_grid_orog"].replace("*GRIDNAME*", gridname) 

img_dpi=int(cfg['plotting']['dpi'])
make_plot=bool(cfg['plotting']['make_plot'])


# Open the file(s) with model gridpoints
lat_model=np.loadtxt(path_grid_in+txt_model_grid_lat_in)
lon_model=np.loadtxt(path_grid_in+txt_model_grid_lon_in)
nlat_model=len(lat_model)
nlon_model=len(lon_model)

# Open the file with 1km-resolution orography smoothed to model gridspacing
data_orog_1km_smooth = xr.open_dataset(path_data_in+netcdf_1km_smooth_orog_in)

data_orog_1km_smooth=wrap360(data_orog_1km_smooth,'longitude')

lat_1km = data_orog_1km_smooth.latitude.values
lon_1km = data_orog_1km_smooth.longitude.values
nlat_1km=len(lat_1km)
nlon_1km=len(lon_1km)

# AVERAGE THE 1KM SMOOTH OROGRAPHY TO THE MODEL GRIDBOXES.

# START FROM THE FIRST LATITUDE --> SOUTH POLE

# GLOBO POINTS ARE THE T POINTS, SO POINTS AT LAT +90 AND -90 REPRESENT
# CIRCULAR AREAS AROUND THE POLES

# THE CODE EXPECT THE DATA ORDERED IN LATITUDE FROM -90 TO +90

mean_orog_on_model_grid=np.zeros((len(lat_model),len(lon_model)))

for latindex in range(nlat_model):
    
    
    # SOUTH POLE POINT
    if lat_model[latindex]==-90:
        
        print("South Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex+1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        northboundary=-90+dr # lat boundary of the polar cap
        
        # collect all the "1km points" falling in the polar cap
        
        polar_cap_points_1km=data_orog_1km_smooth.sel(latitude=slice(-90,northboundary))
        
        mean_orog_in_this_model_gridpoint=float(polar_cap_points_1km.elev.mean())
        
        mean_orog_on_model_grid[latindex,:]=mean_orog_in_this_model_gridpoint # all points at -90 have the same value, at every longitude
        
        
        
    # NORTH POLE POINT    
    if lat_model[latindex]==90:
        
        print("North Pole")
        
        dr=abs(lat_model[latindex]-lat_model[latindex-1]) # lat are supposed to be ordered south --> north (-90 --> +90)
        
        southboundary=90-dr # lat boundary of the polar cap
        
        # collect all the "1km points" falling in the polar cap
        
        polar_cap_points_1km=data_orog_1km_smooth.sel(latitude=slice(southboundary,90))
        
        mean_orog_in_this_model_gridpoint=float(polar_cap_points_1km.elev.mean())
        
        mean_orog_on_model_grid[latindex,:]=mean_orog_in_this_model_gridpoint # all points at +90 have the same value, at every longitude



    # THESE ARE "NORMAL" LAT-LON POINTS (NOT ON THE POLE)  
    if lat_model[latindex]>-90 and lat_model[latindex]<90:
        
        print("Latitude: ",lat_model[latindex])
        
        for lonindex in range(nlon_model):
            
            dlat_north=abs(lat_model[latindex]-lat_model[latindex+1])/2
            dlat_south=abs(lat_model[latindex]-lat_model[latindex-1])/2
            
            
            
            # FIRST LONGITUDE POINT (IN A "0-360" GRID)
            if lonindex==0:
                
                # if i am at the first model grid longitude, the dlon_west is this lon +360
                # minus the last lon in the array. I AM ASSUMING LONGITUDES
                # FROM 0 TO 360
                dlon_west=abs(lon_model[lonindex]+360-lon_model[-1])/2 
                dlon_east=abs(lon_model[lonindex]-lon_model[lonindex+1])/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                
                eastboundary =lon_model[lonindex]+dlon_east
                                
                if lon_model[lonindex]-dlon_west>=0:
                    
                    # the grid cell west border is at lon >=0
                    westboundary =lon_model[lonindex]-dlon_west
                    #print("first cell west boundary: ",westboundary)
                    
                    points_1km_inside_model_gridbox=data_orog_1km_smooth.sel(latitude=slice(southboundary,northboundary),longitude=slice(westboundary,eastboundary)) 
                
                else:
                    
                    # the grid cell west border goes in the western hemisphere                    
                    data_orog_1km_selection_NS=data_orog_1km_smooth.sel(latitude=slice(southboundary,northboundary))
                    westboundary=360+(lon_model[lonindex]-dlon_west) # the parenthesis is a negative number
                    #print("first cell west boundary: ",westboundary)

                    points_1km_inside_model_gridbox=data_orog_1km_selection_NS.where((data_orog_1km_selection_NS.longitude>westboundary) | (data_orog_1km_selection_NS.longitude<eastboundary),drop=True)
            
            
            
            # LAST LONGITUDE POINT (IN A "0-360" GRID)
            if lonindex==(nlon_model-1):
                
                # if i am at the last model grid lon, the dlon_east is this lon minus
                # the first lon in the array +360. I AM ASSUMING LONGITUDES
                # FROM 0 TO 360
                dlon_west=abs(lon_model[lonindex]-lon_model[lonindex-1])/2 
                dlon_east=abs(lon_model[lonindex]-(lon_model[0]+360))/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                
                westboundary =lon_model[lonindex]-dlon_west
                
                if lon_model[lonindex]+dlon_east<=360:
                    
                    # the grid cell east border is at lon <=360
                    eastboundary =lon_model[lonindex]+dlon_east
                    #print("last cell east boundary: ",eastboundary)
                    
                    points_1km_inside_model_gridbox=data_orog_1km_smooth.sel(latitude=slice(southboundary,northboundary),longitude=slice(westboundary,eastboundary)) 
                
                else:
                    
                    # the grid cell east border goes in the eastern hemisphere                    
                    data_orog_1km_selection_NS=data_orog_1km_smooth.sel(latitude=slice(southboundary,northboundary))
                    eastboundary=(lon_model[lonindex]+dlon_east)-360 
                    #print("last cell east boundary: ",eastboundary)

                    points_1km_inside_model_gridbox=data_orog_1km_selection_NS.where((data_orog_1km_selection_NS.longitude>westboundary) | (data_orog_1km_selection_NS.longitude<eastboundary),drop=True)
                
                
                   
            # "CENTRAL" (NOT FIRST, NOT LAST) LONGITUDE POINT (IN A "0-360" GRID) ("EASY" POINTS)
            if lonindex>0 and lonindex<nlon_model-1: 
                
                # i am not at the first or last model grid longitude
            
                dlon_west=abs(lon_model[lonindex]-lon_model[lonindex-1])/2 
                dlon_east=abs(lon_model[lonindex]-lon_model[lonindex+1])/2
                
                northboundary=lat_model[latindex]+dlat_north
                southboundary=lat_model[latindex]-dlat_south
                westboundary =lon_model[lonindex]-dlon_west
                eastboundary =lon_model[lonindex]+dlon_east
                
                points_1km_inside_model_gridbox=data_orog_1km_smooth.sel(latitude=slice(southboundary,northboundary),longitude=slice(westboundary,eastboundary))
                
             
            
            #####################################################################################
            ### calculate the mean on the selected 1km data points inside the model grid cell ###
            #####################################################################################
            mean_orog_in_this_model_gridpoint=float(points_1km_inside_model_gridbox.elev.mean())
            mean_orog_on_model_grid[latindex,lonindex]=mean_orog_in_this_model_gridpoint
 
                

### create data array and save to netcdf

mean_orog_on_model_grid_da=xr.DataArray(mean_orog_on_model_grid, coords=[('latitude', lat_model),('longitude', lon_model)])

mean_orog_on_model_grid_da.to_dataset(name = 'elev').to_netcdf(path_data_out+netcdf_model_grid_orog_out)

## plotting
if make_plot:
    oroplot.orography_plot(mean_orog_on_model_grid,path_img_out+img_model_grid_orog_out,img_dpi)

