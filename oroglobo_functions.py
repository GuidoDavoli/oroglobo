#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:48:14 2023

@author: guidodavoli
"""

import numpy as np
from shapely import Polygon
import LatLon23
import oroglobo_parameters as oropar
import rioxarray as rioxr
import xarray as xr
import os
import time
from haversine import haversine
start = time.time()


# IMPORT PARAMETERS
path_data_in=oropar.paths_in["copernicus_90m"]
tif_copernicus_90m_in=oropar.files_in["tif_copernicus_90m"]


def get_copernicus90m_tiles_list_in_grid_box(grid_box_polygon):
    
    """
    this function receives a shapely polygon representing a gridbox
    and returns the list of filenames of the copernicus 90m DEM tiles
    which intersect the gridbox.
    
    the simplest way to do it is to scan all the globe one deg by one deg
    (copernicus tiles are 1x1 degs.) and check intersections with the 
    grid_box_polygon.
    But this is a waste of computational time (a lot of times I have no 
    intersections); to avoid time waste, I check intersections only on a region
    given by the boundaries of the grid_box_polygon + a frame around it.
    the frame is thickness is "framethickness" (in degrees)
    """
    
    minlon,minlat,maxlon,maxlat=grid_box_polygon.bounds
    framethickness=2
    # I AM ASSUMING A -180/+180 GRID
    
    # floor: nearest integer before
    # ceil:  nearest integer after    
    minlon=int(max(np.floor(minlon)-framethickness,-180))
    minlat=int(max(np.floor(minlat)-framethickness,-90))
    maxlon=int(min(np.ceil(maxlon)+framethickness,180))
    maxlat=int(min(np.ceil(maxlat)+framethickness,90))
    
    tiles_filenames_list=[]

    #for southboundary_tile in range(-90,90): # numbers in [-90:89]
    for southboundary_tile in range(minlat,maxlat+1): # numbers in [minlat:maxlat]
            
        northboundary_tile=southboundary_tile+1
        
        #for westboundary_tile in range(-180,180): # numbers in [-180,179]
        for westboundary_tile in range(minlon,maxlon+1): # numbers in [minlon,maxlon]
        
            eastboundary_tile=westboundary_tile+1
            
            tile_polygon=Polygon([(westboundary_tile,southboundary_tile), (eastboundary_tile,southboundary_tile), (eastboundary_tile,northboundary_tile), (westboundary_tile,northboundary_tile)])
            
            # both the input grid box and the tile grid are on a common -90,90;180,180 grid, so i can check intersection (overlap) with shapely
            # I USE OVERLAPS&cCONTAINS&WITHIN INSTEAD OF INTERSECTS IN ORDER TO AVOID "TRUE" IF THE TWO PLYGONS SIMPLY "TOUCHES" EACH OTHER
            # BUT THERE IS NO A TRUE OVERLAP (THAT IS THERE ARE NO POINTS OF THE TILE REALLY INSIDE THE input GRID BOX).
            # THIS IS NEEDED TO AVOID ERRORS IN THE NEXT CALCULATIONS.
            
            # TRY ALSO: INTERSECT AND !TOUCH
            if grid_box_polygon.overlaps(tile_polygon) or grid_box_polygon.contains(tile_polygon) or grid_box_polygon.within(tile_polygon):
                
                # convert to ""4 hemispheres" format and open the file
                                
                tile_SW_vertex=LatLon23.LatLon(lat=southboundary_tile,lon=westboundary_tile)
                lat_string=tile_SW_vertex.to_string('d% %H')[0]
                lon_string=tile_SW_vertex.to_string('d% %H')[1]
                
                hemisphere_NS=lat_string[-1]
                hemisphere_EW=lon_string[-1]
                sector=hemisphere_NS+hemisphere_EW
                path_dem=path_data_in.replace("*SECTOR*", sector)
                
                lat_string=lat_string.replace(hemisphere_NS, '') # remains only the number
                lon_string=lon_string.replace(hemisphere_EW, '') # remains only the number
                
                lat_tile=int(lat_string)
                lon_tile=int(lon_string)
                
                lat_tile_name=str(lat_tile).zfill(2)
                lon_tile_name=str(lon_tile).zfill(3)
                tif_tile_name=tif_copernicus_90m_in
                tif_tile_name=tif_tile_name.replace("*HEMISPHERE-NS*", hemisphere_NS)
                tif_tile_name=tif_tile_name.replace("*HEMISPHERE-EW*", hemisphere_EW)
                tif_tile_name=tif_tile_name.replace("*LAT2DIGITS*", lat_tile_name)
                tif_tile_name=tif_tile_name.replace("*LON3DIGITS*", lon_tile_name)
                tif_tile_name=path_dem+tif_tile_name
                
                tiles_filenames_list.append(tif_tile_name)
                
    return tiles_filenames_list


def filter_lowpass_2d_van_niekerk(lon,lat,data2d):
    
    # CONTA SE DATA2D HA DIMENSIONI LON,LAT OPPURE LAT,LON? CONTROLLA
    
    ### ADAPTED FROM THE CODE PUBLICHED WITH THE PAPER BY VAN NIEKERK, 2021, QJRMS
    
    lat_rad=np.radians(lat)
    lon_rad=np.radians(lon)
    
    radius=6371000 # meters
    Lx=5000
    filter='Gaussian'
    n=2
    
    Nx = lon.shape[0]
    My = lat.shape[0]
    dx = radius*np.cos(lat_rad[int(My/2.)])*np.abs(lon_rad[1] - lon_rad[0])
    dy = radius*np.abs(lat_rad[1] - lat_rad[0])
    hhat = np.fft.fft2(data2d)
    k = 2*np.pi*np.fft.fftfreq(Nx,dx)
    l = 2*np.pi*np.fft.fftfreq(My,dy)
    kk, ll = np.meshgrid(k,l)
    ktot = (kk**2 + ll**2)**0.5
    kmax = 2*np.pi/Lx
    #hhat[0,0] = 0.0  # cut off "0" frequency: is like to subtract the mean value from the whole image
    if filter == 'Gaussian':
        res = np.exp(- (ktot/kmax)**2)
    if filter == 'Butterworth':
        res = 1 / np.sqrt(1 + (ktot / kmax) ** (2 * n))  # nth order butterworth
    if filter == 'ideal':
        res=np.ones(ktot.shape)
        res[np.abs(kk)>=kmax] = 0.0
        res[np.abs(ll)>=kmax] = 0.0
        res[np.abs(ktot)>=kmax] = 0.0
    if filter == '':
        res=1
    hfilt = np.real(np.fft.ifft2(hhat*res))
    
    # SEE THE ORIGINAL CODE TO COMPUTE HAMP HERE
    # no, computed in other function
    
    return hfilt


def calculate_F1_F2_F3_hamp(lon,lat,data2d,dx,dy):
    
    # CONTA SE DATA2D HA DIMENSIONI LON,LAT OPPURE LAT,LON? CONTROLLA

    # preso da codice van niekerk (adapted)
    
    # dx,dy are in METERS and calculated outside the function 
    # (alternative method wrt van niekerk)
    
    hmin = data2d.min()
    hmax = data2d.max()
    ogwd_hamp = hmax - hmin
    
    #### Commented since dx and dy are now calculated outside and received as parameters
    #radius=6371000 # meters
    #lat_rad=np.radians(lat)
    #lon_rad=np.radians(lon)
    Nx = lon.shape[0]
    My = lat.shape[0]
    #dx = radius*np.cos(lat_rad[int(My/2.)])*np.abs(lon_rad[1] - lon_rad[0])
    #dy = radius*np.abs(lat_rad[1] - lat_rad[0])
    
    
    hhat = np.fft.fft2(data2d)*dx*dy/(4*np.pi**2)
    hhat[0,0] = 0
    k = 2*np.pi*np.fft.fftfreq(Nx,dx)
    l = 2*np.pi*np.fft.fftfreq(My,dy)
    kk, ll = np.meshgrid(k,l)
    ktot = (kk**2 + ll**2)**0.5
    dk = np.abs(k[1] - k[0])
    dl = np.abs(l[1] - l[0])
    area = My*dy * Nx*dx
    ogwd_F1 = np.nansum( ((4*np.pi**2 * kk**2 * np.abs(hhat)**2 * dk * dl / ktot) /area) )
    ogwd_F2 = np.nansum( ((4*np.pi**2 * kk*ll * np.abs(hhat)**2 * dk * dl / ktot) /area) )
    ogwd_F3 = np.nansum( ((4*np.pi**2 * ll**2 * np.abs(hhat)**2 * dk * dl / ktot) /area) )
    #sd = np.nansum(((4*np.pi**2 * np.abs(hhat)**2 * dk * dl ) /area) )
    
    return ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp

def calculate_stddev_anis_orient_slope(data2d,dx,dy):
    
    # CONTA SE DATA2D HA DIMENSIONI LON,LAT OPPURE LAT,LON? CONTROLLA
    
    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    
    # https://www.ecmwf.int/sites/default/files/elibrary/2016/16661-experience-creating-orography-ancillary-files.pdf

    
    IN NUMYPY, THE '0' AXIS IS VERTICAL ('ROWS') GOING FROM UP TO BOTTOM (LIKE A '-Y' CARTESIAN AXIS);
               THE '1' AXIS IS HORIZONTAL ('COLUMNS') GOING FROM LEFT TO RIGHT (LIKE THE 'X' CARTESIAN AXIS).
               EXACTLY HOW IMSHOW PLOTS THE AXIS LABELS BY DEFAULT
               
    dx and dy are the lenght of each pixel in the two dimensions (y=lat, x=lon), in METERS
               
    """
    
    stddev=np.std(data2d)

    grad_y=-np.gradient(data2d, dy, axis=0)   ### SEE EXPLANATION BEFORE
    grad_x=np.gradient(data2d, dx, axis=1)
    """
    plt.figure()
    plt.imshow(grad_x,cmap='seismic',vmin=-100,vmax=100)
    plt.colorbar()

    plt.figure()
    plt.imshow(grad_y,cmap='seismic',vmin=-100,vmax=100)
    plt.colorbar()
    """
    K=0.5* ( np.mean(np.square(grad_x)) + np.mean(np.square(grad_y)) )

    L=0.5* ( np.mean(np.square(grad_x)) - np.mean(np.square(grad_y)) )

    M=np.mean( np.multiply(grad_x,grad_y) )
    
    if M==0 or L==0:
        
        orient_rad=np.nan # it means: undefined
    
    if L>0:
        
        orient_rad= 0.5*np.arctan(M/L)   # È RISP AD ASSE X cartesiano, in radianti
        
    if L<0 and M>0:

        orient_rad= 0.5*np.arctan(M/L) + np.pi/2
        
    if L<0 and M<0:

        orient_rad= 0.5*np.arctan(M/L) - np.pi/2   
        
    orientation=np.degrees(orient_rad)


    Lprime=np.sqrt(L**2+M**2)

    anisotropy=np.sqrt( (K-Lprime)/(K+Lprime) )

    slope = np.sqrt( K + np.sqrt( np.square(L) + np.square(M) ) )
    
    
    return stddev,anisotropy,orientation,slope




def calculate_ogwd_and_tofd_parameters_in_model_grid_box(model_grid_box_polygon,all_tif_data_ds,operational_mean_orog_in_model_grid_box):
    
    """
    this function receives:
        - a shapely polygon describing the model grid box
        - an xarray dataset with the data of the DEM tiles intersecting the model grid box
        - the value of the model operational mean orography inside the model grid box
    and returns:
        - the values of all the ogwd parameters for that model grid box
    that is done in the following steps:
        1) create a xarray dataset with the DEM data, ordered by lat and lon;
        if they are not on the disk,
        it means that there is only sea on that area --> fill with zero.
        2) subtract the model operational mean orography, and obtain the so called
            "subgrid orography"
        .
        .
        .
        .
        
    """
        
    # calculate the subgrid orog, that is the deviation of the high res data from the 
    # model operational mean grid 
    subgrid_elev_all_tif_data_ds=all_tif_data_ds-operational_mean_orog_in_model_grid_box
    """
    # plot for debug
    all_tif_data_ds.elev.plot(cmap='terrain',vmin=0,vmax=4000)
    plt.show()
    subgrid_elev_all_tif_data_ds.elev.plot(cmap='RdBu_r',vmin=-200,vmax=200)
    plt.show()
    """
    # filter out scales > 5km    
    lon=subgrid_elev_all_tif_data_ds.x.values
    lat=subgrid_elev_all_tif_data_ds.y.values
    data2d=subgrid_elev_all_tif_data_ds.elev.data[0]
    subgrid_elev_all_tif_data_5kmfilt_lowpass=filter_lowpass_2d_van_niekerk(lon, lat, data2d)
    
    # make a new xarray dataarray with lat, lon, filtered data
    subgrid_elev_all_tif_data_5kmfilt_lowpass_da=xr.DataArray(
        subgrid_elev_all_tif_data_5kmfilt_lowpass,
        coords=[('latitude', lat),('longitude', lon)])
    """
    # plot for debug
    subgrid_elev_all_tif_data_5kmfilt_lowpass_da.plot(cmap='RdBu_r',vmin=-200,vmax=200)
    plt.show()
    """
    # derive scales <5km as a difference between unfiltered and filtered >5km:
    subgrid_elev_all_tif_data_5kmfilt_highpass=subgrid_elev_all_tif_data_ds.elev.data[0]-subgrid_elev_all_tif_data_5kmfilt_lowpass
    
    # make a new xarray dataarray with lat, lon, high pass data
    subgrid_elev_all_tif_data_5kmfilt_highpass_da=xr.DataArray(
        subgrid_elev_all_tif_data_5kmfilt_highpass,
        coords=[('latitude', lat),('longitude', lon)])
    """
    # plot for debug
    subgrid_elev_all_tif_data_5kmfilt_highpass_da.plot(cmap='RdBu_r',vmin=-20,vmax=20)
    plt.show()
    """
    # cut out only the point really inside the model grid box
    model_grid_box_westboundary = model_grid_box_polygon.bounds[0]
    model_grid_box_southboundary= model_grid_box_polygon.bounds[1]
    model_grid_box_eastboundary = model_grid_box_polygon.bounds[2]
    model_grid_box_northboundary= model_grid_box_polygon.bounds[3]
    #print(model_grid_box_westboundary,model_grid_box_eastboundary,model_grid_box_southboundary,model_grid_box_northboundary)
    subgrid_elev_inside_gridbox_5kmfilt_lowpass_da = subgrid_elev_all_tif_data_5kmfilt_lowpass_da.sel(latitude=slice(model_grid_box_northboundary,model_grid_box_southboundary),longitude=slice(model_grid_box_westboundary,model_grid_box_eastboundary))
    subgrid_elev_inside_gridbox_5kmfilt_highpass_da=subgrid_elev_all_tif_data_5kmfilt_highpass_da.sel(latitude=slice(model_grid_box_northboundary,model_grid_box_southboundary),longitude=slice(model_grid_box_westboundary,model_grid_box_eastboundary))
    
    
    lon_inside_gridbox=subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.longitude.values
    lat_inside_gridbox=subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.latitude.values
    
    # determine the resolution in degrees, then convert to meters with haversine
    # in order to give correct dx and dy to functions that calculates
    # the orographic parameters with gradients
    
    dlon_inside_gridbox=abs(lon_inside_gridbox[1]-lon_inside_gridbox[0]) 
    dlat_inside_gridbox=abs(lat_inside_gridbox[1]-lat_inside_gridbox[0])
    
    central_lon_inside_gridbox=np.median(lon_inside_gridbox)
    central_lat_inside_gridbox=np.median(lat_inside_gridbox)
    
    point2_lon=central_lon_inside_gridbox+dlon_inside_gridbox
    point2_lat=central_lat_inside_gridbox+dlat_inside_gridbox
    
    dem_dlat_meters=haversine((central_lat_inside_gridbox,central_lon_inside_gridbox),(point2_lat,central_lon_inside_gridbox),unit='m')
    dem_dlon_meters=haversine((central_lat_inside_gridbox,central_lon_inside_gridbox),(central_lat_inside_gridbox,point2_lon),unit='m')
    
    """    
    # plot for debug
    subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.plot(cmap='RdBu_r',vmin=-200,vmax=200)
    plt.show()
    subgrid_elev_inside_gridbox_5kmfilt_highpass_da.plot(cmap='RdBu_r',vmin=-20,vmax=20)
    plt.show()
    """
    
    # from that points, calculate the orographic parameters for the grid box
    
    # call specific functions for each group of parameters
   
    ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp=calculate_F1_F2_F3_hamp(
                                        lon_inside_gridbox,
                                        lat_inside_gridbox,
                                        subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.data,
                                        dem_dlon_meters,
                                        dem_dlat_meters
                                        )
    
    ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope=calculate_stddev_anis_orient_slope(
                                                                subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.data,
                                                                dem_dlon_meters,
                                                                dem_dlat_meters
                                                                )
    
    tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope=calculate_stddev_anis_orient_slope(
                                                                subgrid_elev_inside_gridbox_5kmfilt_highpass_da.data,
                                                                dem_dlon_meters,
                                                                dem_dlat_meters
                                                                )
    
    
    # return the values of the parameters
    return ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
           ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
           tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope
           
           
           
def calculate_mean_orog_in_grid_box(km1_grid_box_polygon,all_tif_data_ds):
    
    """
    this function receives:
        - a shapely polygon describing the input grid box
        - an xarray dataset with the data of the DEM tiles intersecting the model grid box 
    and returns:
        - the value of the mean orography inside the input grid box
    that is done in the following steps:
        1) create a xarray dataset with the DEM data, ordered by lat and lon;
        if they are not on the disk,
        it means that there is only sea on that area --> fill with zero.
        2) calculate mean orog
        
    """
        
    all_tif_data=all_tif_data_ds.elev.data[0]
    lon=all_tif_data_ds.x.values
    lat=all_tif_data_ds.y.values
    # make a new xarray dataarray with lat, lon, orog
    all_tif_data_da=xr.DataArray(
       all_tif_data,
       coords=[('latitude', lat),('longitude', lon)])
    
    # cut out only the point really inside the 1km grid box
    km1_grid_box_westboundary = km1_grid_box_polygon.bounds[0]
    km1_grid_box_southboundary= km1_grid_box_polygon.bounds[1]
    km1_grid_box_eastboundary = km1_grid_box_polygon.bounds[2]
    km1_grid_box_northboundary= km1_grid_box_polygon.bounds[3]
    #print(model_grid_box_westboundary,model_grid_box_eastboundary,model_grid_box_southboundary,model_grid_box_northboundary)
    all_tif_data_inside_1km_gridbox_da = all_tif_data_da.sel(latitude=slice(km1_grid_box_northboundary,km1_grid_box_southboundary),longitude=slice(km1_grid_box_westboundary,km1_grid_box_eastboundary))
    
    meanorog=np.mean(all_tif_data_inside_1km_gridbox_da.data)
        
    return meanorog



def get_copernicus90m_data_xrds_from_tiles_list(tif_filenames_list):
    
    """
    This functions receives:
        - a list of filenames corresponding to tiles of the copernicus 90m dataset
    and returns:
        - an xarray dataset containing the copernicus90m data contained
        in the list of tiles filenames, if data exist, and data_exists=1
        - data_exist=0 and an empty list if data are not present
    """
    
    data_exist=1
    
    tif_data_list=[]
    tif_dlon_list=[]
    
    for tif_filename in tif_filenames_list:
        
        if os.path.isfile(tif_filename):
            tif_data = rioxr.open_rasterio(tif_filename)[:,:-1,:-1].to_dataset(name='elev')           
            # if the tile is not on disk, it means that it is all on sea --> do not consider it
            # in the following, missing tiles will be treated as elev=0 
            tif_dlon=abs(tif_data.x.values[1]-tif_data.x.values[0]) # this can vary across the copernicus dataset
            tif_dlat=abs(tif_data.y.values[1]-tif_data.y.values[0]) # this is always the same, data are uniform in latitude
            
            tif_data_list.append(tif_data)
            tif_dlon_list.append(tif_dlon)

    if len(tif_data_list)>0: # if at least one tile exist --> if there is land, i have data
    
        # i have to check if the selected tif data have all the same "dlon"
        # (i.e. the same longitudinal coordinates). Copernicus data have different
        # dlon in different latitudinal bands in order to mantain an approx. costant
        # spacing in meters. So when the model gridbox crosses specific latitudes
        # (like +-85, +-80, +-70...+-50) tif files with different longitudinal
        # coordinates and spacing are loaded; they must be "homogenized"
        # before passing it to xr.combine_by_coords
        
        # are all dlons equal or not?
        
        if ( tif_dlon_list.count(tif_dlon_list[0]) == len(tif_dlon_list) ) == False:
            
            # if they are not all equal, proceed by steps:
                
            # 1) combine with combine_by_coords and fill with nan where there are no data. 
            #    this means that will became nan:
            #    - points on the sea
            #    - the "spurious/duplicate" points that arise beacuse the grids are not equal 
            #      (combine_by_coords create a grid with all longitudes in the two grids)  
            all_tif_data_ds=xr.combine_by_coords(tif_data_list,fill_value=np.nan)
            
            # 2) fill the nans forward and backward with the nearest value after and before (limit="1 pixel")
            #   (since in the worst case the fine grid spacing is half of the coarse grid spacing, it is supposed that limit=1 before and after is enough to fill gaps  
            #   some sea points (the closest to the coasts) are "lost" (i.e. the became >0) with this procedure 
            all_tif_data_ds=all_tif_data_ds.bfill(dim='x',limit=1)
            all_tif_data_ds=all_tif_data_ds.ffill(dim='x',limit=1)
            
            # 3) transform all the nan still present to zero (i suppose that the remaining nans are sea)
            all_tif_data_ds=all_tif_data_ds.fillna(0)
            
            # 4) get the larger dlon (lower resolution)
            tif_larger_dlon=np.array(tif_dlon_list).max()
            
            # 5) get the bounds of the tiff data all together
            all_tiff_westboundary= all_tif_data_ds.x.values[0]
            all_tiff_eastboundary= all_tif_data_ds.x.values[-1]
            all_tiff_northboundary=all_tif_data_ds.y.values[0]
            all_tiff_southboundary=all_tif_data_ds.y.values[-1]
            
            # 6) build lat and lon arrays, equally spaced with spacing tif_larger_dlon,
            #    which will be used to resample the data
            newlat_resample=np.arange(all_tiff_northboundary,all_tiff_southboundary-0.001,-tif_dlat) # +0.001 to not exclude the last point
            newlon_resample=np.arange(all_tiff_westboundary,all_tiff_eastboundary+0.001,tif_larger_dlon) # +0.001 to not exclude the last point
            
            # 7) resample the data at uniform resolution, consistent with the 
            # lowest resolution tile
            
            all_tif_data_ds=all_tif_data_ds.interp(coords={"x": newlon_resample, "y": newlat_resample})
            # again set nan to zero, if any
            all_tif_data_ds=all_tif_data_ds.fillna(0)
            
        else:
            
            # if yes, data are have uniform grids, so
            # exploit xarray to order the tiles by lat and lon; simply fill with zero where there are no data 
            # it happens olny on the sea)
            all_tif_data_ds=xr.combine_by_coords(tif_data_list,fill_value=0) # where no data (sea) --> =0
        
    else: # if there are no data on disk --> the grid box is entirely on ocean
        
        all_tif_data_ds=[]
        data_exist=0
    
    return all_tif_data_ds,data_exist