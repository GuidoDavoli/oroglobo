#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:48:14 2023

@author: guidodavoli
"""

import numpy as np
import math
from shapely import Polygon
import LatLon23
import oroglobo_parameters as oropar
import rioxarray as rioxr
import xarray as xr
import os
import matplotlib.pyplot as plt


# IMPORT PARAMETERS
path_data_in=oropar.paths_in["copernicus_90m"]
tif_copernicus_90m_in=oropar.files_in["tif_copernicus_90m"]


def distance(a,b):
    
    dist=np.sqrt(a**2 + b**2)
    
    return dist


def smoothing2D_ECMWF(r, d, D):
    
    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
    R = np.absolute(r)
    
    a = D/2 - d
    b = D/2 + d
    
    if R < a:
        
        h = 1/D
        
    if a <= R and R <= b:
        
        h = 1/(2*D) + 1/(2*D)*np.cos( np.pi*(r-D/2+d) )/(2*d)
        
    if R > b:
        
        h = 0
        
    return h


def filter_ECMWF(d,D):


    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    """
    
    #### prepare the 2d square array representing the filter (scale in pixel)
    #### the filter is a square with an odd number of pixels (L)
    #### but L is determinig in different ways depending on D (even or odd)
    
    if (D%2)==0: # D is even
    
        L=int((D/2+d)*2+1)
        
    else:
        
        L=int(math.ceil(D/2+d)*2+1) # ceil(x) = round to the smallest integer higher or equal to x 
        
    
    filt = np.zeros((L,L))
    
    i0_filt=-int((L-1)/2)
    j0_filt=-int((L-1)/2)
    
    for i in range(L):
        for j in range(L):
            
            dist=distance(i0_filt+i, j0_filt+j)
            filt[i,j]=smoothing2D_ECMWF(dist, d, D)
            
    filt=filt/np.sum(filt)  #### tentativo di conservare il "modulo" dell'orografia
                            #### ispirato da https://medium.com/@bdhuma/6-basic-things-to-know-about-convolution-daef5e1bc411
                            #### the idea is to put a normalization factor in front of the filter matrix
    
    return filt


def get_copernicus90m_tiles_list_in_model_grid_box(model_grid_box_polygon):
    
    """
    this function receives a shapely polygon representing the model gridbox
    and returns the list of filenames of the copernicus 90m DEM tiles
    which intersect the model gridbox
    """
    
    tiles_filenames_list=[]

    for southboundary_tile in range(-90,90): # numbers in [-90:89]
            
        northboundary_tile=southboundary_tile+1

        for westboundary_tile in range(-180,180): # numbers in [-180,179]
        
            eastboundary_tile=westboundary_tile+1
            
            tile_polygon=Polygon([(westboundary_tile,southboundary_tile), (eastboundary_tile,southboundary_tile), (eastboundary_tile,northboundary_tile), (westboundary_tile,northboundary_tile)])
            
            # both the model grid and the tile grid are on a common -90,90;180,180 grid, so i can check intersection with shapely
            if model_grid_box_polygon.intersects(tile_polygon):
                
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
    
    return hfilt


def calculate_F1_F2_F3_hamp():
    
    # prendi da codice van niekerk
    
    return

def calculate_stddev_anis_orient_slope(data2d):
    
    """
    SEE https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
    
    # https://www.ecmwf.int/sites/default/files/elibrary/2016/16661-experience-creating-orography-ancillary-files.pdf

    
    IN NUMYPY, THE '0' AXIS IS VERTICAL ('ROWS') GOING FROM UP TO BOTTOM (LIKE A '-Y' CARTESIAN AXIS);
               THE '1' AXIS IS HORIZONTAL ('COLUMNS') GOING FROM LEFT TO RIGHT (LIKE THE 'X' CARTESIAN AXIS).
               EXACTLY HOW IMSHOW PLOTS THE AXIS LABELS BY DEFAULT
    """
    
    stddev=np.std(data2d)

    grad_y=-np.gradient(data2d, axis=0)   ### SEE EXPLANATION BEFORE
    grad_x=np.gradient(data2d, axis=1)

    plt.figure()
    plt.imshow(grad_x,cmap='seismic',vmin=-100,vmax=100)
    plt.colorbar()

    plt.figure()
    plt.imshow(grad_y,cmap='seismic',vmin=-100,vmax=100)
    plt.colorbar()

    K=0.5* ( np.mean(np.square(grad_x)) + np.mean(np.square(grad_y)) )

    L=0.5* ( np.mean(np.square(grad_x)) - np.mean(np.square(grad_y)) )

    M=np.mean( np.multiply(grad_x,grad_y) )

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



def calculate_ogwd_and_tofd_parameters_in_model_grid_box(model_grid_box_polygon,tif_filenames_list,operational_mean_orog_in_model_grid_box):
    
    """
    this function receives:
        - a shapely polygon describing the model grid box
        - a list of the DEM tiles intersecting the model grid box
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
    
    tif_data_list=[]
    
    for tif_filename in tif_filenames_list:
        
        if os.path.isfile(tif_filename):
            tif_data = rioxr.open_rasterio(tif_filename)[:,:-1,:-1].to_dataset(name='elev')           
            # if the tile is not on disk, it means that it is all on sea --> do not consider it
            # in the following, missing tiles will be treated as elev=0 
            tif_data_list.append(tif_data)

    # exploit xarray to order the tiles by lat and lon 
    all_tif_data_ds=xr.combine_by_coords(tif_data_list,fill_value=0) # where no data (sea) --> =0
    
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
    """    
    # plot for debug
    subgrid_elev_inside_gridbox_5kmfilt_lowpass_da.plot(cmap='RdBu_r',vmin=-200,vmax=200)
    plt.show()
    subgrid_elev_inside_gridbox_5kmfilt_highpass_da.plot(cmap='RdBu_r',vmin=-20,vmax=20)
    plt.show()
    """
    
    # from that points, calculate the orographic parameters for the grid box
    
        # call specific functions for each parameter
    
    X=(model_grid_box_westboundary+
        model_grid_box_eastboundary+
        model_grid_box_southboundary+
        model_grid_box_northboundary)/4 # just a crazy number to assign as a test
   
    ogwd_F1=X
    ogwd_F2=X
    ogwd_F3=X
    ogwd_hamp=X
    
    ogwd_stddev=X
    ogwd_anisotropy=X
    ogwd_orientation=X
    ogwd_slope=X
    
    tofd_stddev=X
    tofd_anisotropy=X
    tofd_orientation=X
    tofd_slope=X
    
    # return the values of the parameters
    
    return ogwd_F1,ogwd_F2,ogwd_F3,ogwd_hamp,\
           ogwd_stddev,ogwd_anisotropy,ogwd_orientation,ogwd_slope,\
           tofd_stddev,tofd_anisotropy,tofd_orientation,tofd_slope