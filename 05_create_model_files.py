#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

This code:

    - perform the last operations needed to store OroGlobo output
    in netcdf files with a structure fully compatible with the GLOBO model.

"""

import numpy as np
import xarray as xr
import oroglobo_parameters as oropar
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


######################## IMPORT PARAMETERS ##########################

configname='oroglobo_parameters.yaml'
with open(configname, 'r', encoding='utf-8') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

gridname=cfg['model_grid']['GRIDNAME']

path_data_out=cfg['paths_out']["data_out"].replace("*GRIDNAME*", gridname) 


############# open and read tofd file

data_ds=xr.open_dataset(path_data_out+'model_grid_tofd_params_'+gridname+'.nc')

data_ds=wrap360(data_ds,'longitude')

tofd_stddev=data_ds.tofd_stddev.data

new_tofd_stddev=np.zeros((tofd_stddev.shape[0],tofd_stddev.shape[1]+2))

new_tofd_stddev[:,1:-1]=tofd_stddev # fill with stddev, except ghost points

new_tofd_stddev[:,0]=new_tofd_stddev[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_tofd_stddev[:,-1]=new_tofd_stddev[:,1] # repeat extra longitudes. like in geo.F90 row 465

new_tofd_stddev=np.nan_to_num(new_tofd_stddev)  # replace nan with zeros


############# open and read ogwd file

data_ds=xr.open_dataset(path_data_out+'model_grid_ogwd_params_'+gridname+'.nc')

data_ds=wrap360(data_ds,'longitude')

ogwd_stddev=data_ds.ogwd_stddev.data
ogwd_slope=data_ds.ogwd_slope.data
ogwd_orientation=data_ds.ogwd_orientation.data
ogwd_anisotropy=data_ds.ogwd_anisotropy.data
ogwd_hamp=data_ds.ogwd_hamp.data
ogwd_F1=data_ds.ogwd_F1.data
ogwd_F2=data_ds.ogwd_F2.data
ogwd_F3=data_ds.ogwd_F3.data

new_ogwd_stddev=np.zeros((ogwd_stddev.shape[0],ogwd_stddev.shape[1]+2))
new_ogwd_slope=np.zeros((ogwd_slope.shape[0],ogwd_slope.shape[1]+2))
new_ogwd_orientation=np.zeros((ogwd_orientation.shape[0],ogwd_orientation.shape[1]+2))
new_ogwd_anisotropy=np.zeros((ogwd_anisotropy.shape[0],ogwd_anisotropy.shape[1]+2))
new_ogwd_hamp=np.zeros((ogwd_hamp.shape[0],ogwd_hamp.shape[1]+2))
new_ogwd_F1=np.zeros((ogwd_F1.shape[0],ogwd_F1.shape[1]+2))
new_ogwd_F2=np.zeros((ogwd_F2.shape[0],ogwd_F2.shape[1]+2))
new_ogwd_F3=np.zeros((ogwd_F3.shape[0],ogwd_F3.shape[1]+2))

new_ogwd_stddev[:,1:-1]=ogwd_stddev # fill with stddev, except ghost points
new_ogwd_slope[:,1:-1]=ogwd_slope # fill with         , except ghost points
new_ogwd_orientation[:,1:-1]=ogwd_orientation # fill with         , except ghost points
new_ogwd_anisotropy[:,1:-1]=ogwd_anisotropy # fill with         , except ghost points
new_ogwd_hamp[:,1:-1]=ogwd_hamp # fill with         , except ghost points
new_ogwd_F1[:,1:-1]=ogwd_F1 # fill with         , except ghost points
new_ogwd_F2[:,1:-1]=ogwd_F2 # fill with         , except ghost points
new_ogwd_F3[:,1:-1]=ogwd_F3 # fill with         , except ghost points

new_ogwd_stddev[:,0]=new_ogwd_stddev[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_stddev[:,-1]=new_ogwd_stddev[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_slope[:,0]=new_ogwd_slope[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_slope[:,-1]=new_ogwd_slope[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_orientation[:,0]=new_ogwd_orientation[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_orientation[:,-1]=new_ogwd_orientation[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_anisotropy[:,0]=new_ogwd_anisotropy[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_anisotropy[:,-1]=new_ogwd_anisotropy[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_hamp[:,0]=new_ogwd_hamp[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_hamp[:,-1]=new_ogwd_hamp[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F1[:,0]=new_ogwd_F1[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F1[:,-1]=new_ogwd_F1[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F2[:,0]=new_ogwd_F2[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F2[:,-1]=new_ogwd_F2[:,1] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F3[:,0]=new_ogwd_F3[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_ogwd_F3[:,-1]=new_ogwd_F3[:,1] # repeat extra longitudes. like in geo.F90 row 465

new_ogwd_stddev=np.nan_to_num(new_ogwd_stddev)  # replace nan with zeros
new_ogwd_slope=np.nan_to_num(new_ogwd_slope)  # replace nan with zeros
new_ogwd_orientation=np.nan_to_num(new_ogwd_orientation)  # replace nan with zeros
new_ogwd_anisotropy=np.nan_to_num(new_ogwd_anisotropy)  # replace nan with zeros
new_ogwd_hamp=np.nan_to_num(new_ogwd_hamp)  # replace nan with zeros
new_ogwd_F1=np.nan_to_num(new_ogwd_F1)  # replace nan with zeros
new_ogwd_F2=np.nan_to_num(new_ogwd_F2)  # replace nan with zeros
new_ogwd_F3=np.nan_to_num(new_ogwd_F3)  # replace nan with zeros


############# open and read mean orography file

data_ds=xr.open_dataset(path_data_out+'model_grid_operational_orog_'+gridname+'.nc')

data_ds=wrap360(data_ds,'longitude')

elev=data_ds.elev.data

new_elev=np.zeros((elev.shape[0],elev.shape[1]+2))

new_elev[:,1:-1]=elev # fill with elev, except ghost points

new_elev[:,0]=new_elev[:,-2] # repeat extra longitudes. like in geo.F90 row 465
new_elev[:,-1]=new_elev[:,1] # repeat extra longitudes. like in geo.F90 row 465

new_elev=np.nan_to_num(new_elev)  # replace nan with zeros

############ save netcdf

lats=data_ds.latitude.data

dlon=data_ds.longitude.data[1]
lons=np.zeros(len(data_ds.longitude.data)+2)
lons[:-2]=data_ds.longitude.data
lons[-2]=data_ds.longitude.data[-1]+dlon
lons[-1]=data_ds.longitude.data[-1]+dlon+dlon

new_tofd_stddev_da=xr.DataArray(new_tofd_stddev, coords=[('latitude', lats),('longitude', lons)], name="sgo_stddev_tofd")

new_ogwd_stddev_da=xr.DataArray(new_ogwd_stddev, coords=[('latitude', lats),('longitude', lons)], name="sgo_stddev_ogwd")
new_ogwd_slope_da=xr.DataArray(new_ogwd_slope, coords=[('latitude', lats),('longitude', lons)], name="sgo_slope_ogwd")
new_ogwd_orientation_da=xr.DataArray(new_ogwd_orientation, coords=[('latitude', lats),('longitude', lons)], name="sgo_orientation_ogwd")
new_ogwd_anisotropy_da=xr.DataArray(new_ogwd_anisotropy, coords=[('latitude', lats),('longitude', lons)], name="sgo_anisotropy_ogwd")
new_ogwd_hamp_da=xr.DataArray(new_ogwd_hamp, coords=[('latitude', lats),('longitude', lons)], name="sgo_hamp_ogwd")
new_ogwd_F1_da=xr.DataArray(new_ogwd_F1, coords=[('latitude', lats),('longitude', lons)], name="sgo_F1_ogwd")
new_ogwd_F2_da=xr.DataArray(new_ogwd_F2, coords=[('latitude', lats),('longitude', lons)], name="sgo_F2_ogwd")
new_ogwd_F3_da=xr.DataArray(new_ogwd_F3, coords=[('latitude', lats),('longitude', lons)], name="sgo_F3_ogwd")

new_elev_da=xr.DataArray(new_elev, coords=[('latitude', lats),('longitude', lons)], name="elev")


sgo_ds=xr.merge((new_tofd_stddev_da,
               new_ogwd_stddev_da,
               new_ogwd_slope_da,
               new_ogwd_orientation_da,
               new_ogwd_anisotropy_da,
               new_ogwd_hamp_da,
               new_ogwd_F1_da,
               new_ogwd_F2_da,
               new_ogwd_F3_da
               ))

elev_ds=new_elev_da.to_dataset(name = 'elev')

sgo_ds.to_netcdf('oroglobo_rsync/'+'sgo_'+gridname+'_oroglobo.nc',engine="h5netcdf")
elev_ds.to_netcdf('oroglobo_rsync/'+'elev_'+gridname+'_oroglobo.nc',engine="h5netcdf")


