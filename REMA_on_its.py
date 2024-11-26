# -*- coding: utf-8 -*-
"""
Code to interpolate the REMA data onto the ITSLIVE grid
"""

import numpy as np
#import itertools
from netCDF4 import Dataset
from scipy.interpolate import griddata 
#from scipy.fft import fft2, ifft2
#from file_processing import Write_to_nc
#from inversion_module_v2 import smooth_data_load, filter, calcLM_KK, smooth_arr, bed_conditions_clean
#import matplotlib.pyplot as plt
#import time
#from shapely.geometry import box
#import geopandas as gpd

def SURF_coord_NC(map_X, map_Y, map_B, filename):
    ncfile = Dataset(filename,mode='w',format='NETCDF4')
    ncfile.title = 'Output Bed Elevation'
    # Make data dimensions
    ncfile.createDimension('x', map_X.shape[1])     # x axis
    ncfile.createDimension('y', map_X.shape[0])     # y axis
    # Create the axes/variable lengths
    x = ncfile.createVariable('x', np.float32, ('x',))
    x.units = 'km (stereographic) EW'
    y = ncfile.createVariable('y', np.float32, ('y',))
    y.units = 'km (stereographic) NS'
    # Define 2D variables to hold the data
    bed = ncfile.createVariable('Band1',np.float64,('y','x'))
    bed.units = 'm'
    bed.standard_name = 'Inversion output bed'
    # Write in the data
    x[:] = map_X[0,:]
    y[:] = map_Y[:,0]
    bed[:,:] = map_B[:,:].real
    ncfile.close()
    
print('Welcome to this interpolation code')

## Select the data files to be processed
filepath_rema = '../Data/GaplessREMA100.nc'
filepath_itslive = "../Data/antarctica_ice_velocity_450m_v2.nc"

bounds = [-2600000, 2700000, -2100000, 2200000]
# For the rema data
fh = Dataset(filepath_rema, 'r', format='NETCDF4');
X = fh.variables['x'][:]
Y = fh.variables['y'][:]
xl = next(x for x, val in enumerate(X) if val >= bounds[0]) # X in ascending order
xh = next(x for x, val in enumerate(X) if val >= bounds[1])
yl = next(x for x, val in enumerate(Y) if val >= bounds[2])  # Y in descending order
yh = next(x for x, val in enumerate(Y) if val >= bounds[3])  
SURF = fh.variables['Band1'][yl:yh, xl:xh]
X_rema = fh.variables['x'][xl:xh]
Y_rema = fh.variables['y'][yl:yh]
#X_rema, Y_rema = np.meshgrid(X_rema, Y_rema)
fh.close()
fh = Dataset(filepath_itslive, 'r', format='NETCDF4');
X = fh.variables['x'][:]
Y = fh.variables['y'][:]
xl = next(x for x, val in enumerate(X) if val >= bounds[0]) # X in ascending order
xh = next(x for x, val in enumerate(X) if val >= bounds[1])
yl = next(x for x, val in enumerate(Y) if val <= bounds[2])  # Y in descending order
yh = next(x for x, val in enumerate(Y) if val <= bounds[3])  
VX = fh.variables['VX'][yh:yl, xl:xh]
VY = fh.variables['VY'][yh:yl, xl:xh]
X_its = fh.variables['x'][xl:xh]
Y_its = fh.variables['y'][yh:yl]
X_its, Y_its = np.meshgrid(X_its, Y_its)
fh.close()

print('Data loaded')
print(np.min(X_rema), np.max(X_rema), np.min(Y_rema), np.max(Y_rema))
if np.min(X_its) > np.min(X_rema):
    if np.max(X_its) < np.max(X_rema):
        if np.min(Y_its) > np.min(Y_rema):
            if np.max(Y_its) < np.max(Y_rema):
                print('Yes')

X_rema_flat = np.ndarray.flatten(X_rema)
Y_rema_flat = np.ndarray.flatten(Y_rema)
SURF_flat = np.ndarray.flatten(SURF)
print(np.sum(np.isnan(SURF)))
print(len(X_rema_flat), len(Y_rema_flat), len(SURF_flat), X_its.shape)

print(len(X_rema)*len(Y_rema), SURF.shape[0], SURF.shape[1], SURF.shape[0] * SURF.shape[1])

from scipy.interpolate import interpn
points = (Y_rema, X_rema)
xi = (np.flip(Y_its, axis = 0), X_its)
xi = (Y_its, X_its)
SURF_interp = interpn(points, SURF, xi)#, method = 'cubic'), \
 #                          bounds_error = False, fill_value = np.nan)   
print('Interpolation complete')
print(SURF_interp.shape)
print(np.sum(np.isnan(SURF_interp)))

print(Y_its), print(np.flip(Y_its, axis = 0))

# Save the results
#SURF_coord_NC(X_its, np.flip(Y_its, axis = 0), SURF_interp, '../Data/GaplessREMA100_its.nc')
SURF_coord_NC(X_its, Y_its, SURF_interp, '../Data/GaplessREMA100_its2.nc')
print('Files saved')

# Try to import part of this new file to see how it goes 
bounds = [1375100.0, 1424600.0, -774800.0, -725300.0]
fh = Dataset('../Data/GaplessREMA100_its.nc', 'r', format='NETCDF4');
X = fh.variables['x'][:]
Y = fh.variables['y'][:]
xl = next(x for x, val in enumerate(X) if val >= bounds[0]) # X in ascending order
xh = next(x for x, val in enumerate(X) if val >= bounds[1])
yl = next(x for x, val in enumerate(Y) if val >= bounds[2])  # Y in descending order
yh = next(x for x, val in enumerate(Y) if val >= bounds[3])  
SURF = fh.variables['Band1'][yl:yh, xl:xh]
X_rema = fh.variables['x'][xl:xh]
Y_rema = fh.variables['y'][yl:yh]
X_rema, Y_rema = np.meshgrid(X_rema, Y_rema)
fh.close()
print(X,Y, SURF.shape)