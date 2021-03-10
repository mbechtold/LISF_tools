
import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt
from datetime import datetime, timedelta
import logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
from collections import OrderedDict
from netCDF4 import Dataset, date2num
import copy
from validation_good_practice.ancillary import metrics
import matplotlib.pyplot as plt

def lonlat2colrow(lons,lats, clon, clat):
    """ Find nearest tile (col/row) from any given lon/lat """
    londif = np.abs(lons - clon)
    latdif = np.abs(lats - clat)
    col = np.where(np.abs(londif-londif.min())<0.0001)[0][0]
    row = np.where(np.abs(latdif-latdif.min())<0.0001)[0][0]

    return col, row

def create_nc_cube(root="",parameterfile='lis_input.nc',outtype="ROUTING",date_from="",date_to="",timeres="daily"):

    """ create 3D LIS cube for
    (1) surface variables
    (2) routing variables
    """
    date_from = datetime.strptime(date_from,'%Y-%m-%d')
    date_to = datetime.strptime(date_to,'%Y-%m-%d')
    numdays = (date_to - date_from)
    numdays = numdays.days
    date_list = [date_from + dt.timedelta(days=x) for x in range(numdays)]
    cfullpath = os.path.join(root,outtype)

    cdate=date_list[1]
    cfolder = cdate.strftime("%Y%m")
    ds_ = xr.open_dataset(cfullpath+'/'+cfolder+'/LIS_HIST_'+cdate.strftime("%Y%m%d")+'0000.d01.nc')
    ds = xr.open_dataset(parameterfile)
    lons = ds.variables['lon'].values[0,:]
    lons = np.round(np.unique(lons)[np.isnan(np.unique(lons))==False],4)
    lats = ds.variables['lat'].values[:,0]
    lats = np.round(np.unique(lats)[np.isnan(np.unique(lats))==False],4)
    dimensions = OrderedDict([('time',date_list),('lat',lats), ('lon',lons)])
    fname_out = root+'/'+outtype+'.nc'
    ds = ncfile_init(fname_out, dimensions, list(ds_.keys())[2:])

    for c,cdate in enumerate(date_list):
        cfolder = cdate.strftime("%Y%m")
        ds_ = xr.open_dataset(cfullpath+'/'+cfolder+'/LIS_HIST_'+cdate.strftime("%Y%m%d")+'0000.d01.nc')
        for cvar in list(ds_.variables.keys())[3:]:
            # just surface layer, change to make ncfile_init more flexible.
            if cvar=='SoilMoist_inst' or cvar=='SoilMoist_tavg' or cvar=='SoilTemp_inst' or cvar=='SoilTemp_tavg':
                ds.variables[cvar][c,:,:] = ds_.variables[cvar][0,:,:]
            else:
                ds.variables[cvar][c,:,:] = ds_.variables[cvar]

    ds.close()

def ncfile_init(fname, dimensions, variables):
    """ initiale nc file for 3D cube from multiple lis files """
    ds = Dataset(fname, mode='w')
    timeunit = 'hours since 2000-01-01 00:00'

    # Initialize dimensions
    chunksizes = []
    for dim in dimensions:

        # convert pandas Datetime Index to netCDF-understandable numeric format
        if dim == 'time':
            try:
                dimensions[dim] = date2num(dimensions[dim].to_pydatetime(), timeunit).astype('int32')
            except:
                dimensions[dim] = date2num(pd.to_datetime(dimensions[dim]).to_pydatetime(),timeunit).astype('int32')

        # Files are per default image chunked
        if dim in ['lon','lat']:
            chunksize = len(dimensions[dim])
        else:
            chunksize = 1
        chunksizes.append(chunksize)

        dtype = dimensions[dim].dtype
        ds.createDimension(dim, len(dimensions[dim]))
        ds.createVariable(dim,dtype,
                          dimensions=(dim,),
                          chunksizes=(chunksize,),
                          zlib=True)
        ds.variables[dim][:] = dimensions[dim]

    # Coordinate attributes following CF-conventions
    if "time" in ds.variables:
        ds.variables['time'].setncatts({'long_name': 'time',
                                        'units': timeunit})
    ds.variables['lon'].setncatts({'long_name': 'longitude',
                                   'units':'degrees_east'})
    ds.variables['lat'].setncatts({'long_name': 'latitude',
                                    'units':'degrees_north'})

    # Initialize variables
    for var in variables:
        print("create "+var)
        ds.createVariable(var, 'float32',
                          dimensions=list(dimensions.keys()),
                          chunksizes=chunksizes,
                          fill_value=-9999.,
                          zlib=True)

    return ds

