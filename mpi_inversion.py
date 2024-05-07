## Import packages required for running code
from mpi4py import MPI
import numpy as np
import itertools
from netCDF4 import Dataset
#import geopandas as gpd
from scipy.interpolate import griddata 
from scipy.fft import fft, ifft, fft2, ifft2, fftshift, fftfreq
from file_processing import Write_to_nc

from inversion_module_v2 import smooth_arr, calcLM_KK, filter, tapering_func, \
            smooth_data_load, bed_conditions_clean, terminal_inversion_smooth

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

## Select the data files to be processed

filepath_rema = '../Data/GaplessREMA100.nc'
filepath_bedmach = "../Data/BedMachineAntarctica_2019-11-05_v01.nc"
filepath_itslive = "../Data/antarctica_ice_velocity_450m_v2.nc"

## Select the centre coordinate
#centre_coord = [-1.586e6,-0.122e6]                               # PIG hill istar 7
#centre_coord = [-1542249.2233862362, -214786.63862677815]        # PIG 17
#centre_coord = [-1551063.1745184413, -248155.27206390776]        # PIG 18
#centre_coord = [-1.41e6, -0.455e6]                               # Lower thwaites
#centre_coord = [-1.08e6, -0.455e6]                                # WAIS
#centre_coord = [-1.29e6, -0.455e6]                               # Upperish thwaites
#centre_coord = [-1767352.284231, 187790.347215]                  # Ferrigno
#centre_coord = [-1198530, -20920]                                # Lake Elsworth

## Region to invert over
square_size = 50000                       # In metres
interp_grid_spacing = 100                 # In metres (could be set in code by data if wanted)
tapering = 0.1                            # Remove outer XX% due to tapering

## Physical system properties
m = 1                                     # Sliding law constant
C = 100                                   # Mean slipperiness

## Tuneable parameters
p_b = -2                                  # Smooth filter on wavelengths
p_c = -2                                  # Smooth filter on wavelengths
erB = 1 * (10 **(-3))                     # Least squares balance between non-dimensionalised velocity and surface data
erC = 1 * (10 **(-3))                     # Least squares balance between non-dimensionalised velocity and surface data
CutB = 10                                 # Angle to flow of features removed in inversion results
CutC = 15                                 # Angle to flow of features removed in inversion results
wavcutB = 1                               # Factor of h_bar, wavelength of features removed in inversion results
wavcutC = 2                               # Factor of h_bar, wavelength of features removed in inversion results

## How big is the grid 
n = 3                                     # Number of overlapping grids
adj = [6,6]                               # Number of adjacent grids
centre_include = 1                        # How much of the central grid is included (this feature is now not necessary)
filename = 'Output/clean_code_test_wais.nc'    # Filename of output

trans_funcs = 2003                        # Which version of the transfer functions to use (2003 or 2008)


# Import the list of central coordinates 
fh = Dataset('EAnt_coords_list_x_redo.nc', format = 'NETCDF4')
coords_x = fh.variables['data'][:]
fh.close()
fh = Dataset('EAnt_coords_list_y_redo.nc', format = 'NETCDF4')
coords_y = fh.variables['data'][:]
fh.close()

# Reshape
centre_coords = ([])
for i in range(len(coords_x)):
#for i in range(4):
    centre_coords.append([coords_x[i], coords_y[i]])
print(len(centre_coords), 'coordinate points')
    
if rank == 0:
    data = centre_coords

    # determine the size of each sub-task
    ave, res = divmod(len(data), nprocs)
    counts = [ave + 1 if p < res else ave for p in range(nprocs)]

    # determine the starting and ending indices of each sub-task
    starts = [sum(counts[:p]) for p in range(nprocs)]
    ends = [sum(counts[:p+1]) for p in range(nprocs)]

    # converts data into a list of arrays 
    data = [data[starts[p]:ends[p]] for p in range(nprocs)]
else:
    data = None

data = comm.scatter(data, root=0)

filename_base = 'EA_grids/EA_runs_'
#for i in range(len(coords)):

#for rank in nprocs:
for i in range(len(data)):
    if data[i][0]/1000 < 0:
        x_coord = 'm'+str(abs(int(data[i][0]/1000)))
    else: 
        x_coord = str(abs(int(data[i][0]/1000)))
    if data[i][1]/1000 < 0:
        y_coord = 'm'+str(abs(int(data[i][1]/1000)))
    else: 
        y_coord = str(abs(int(data[i][1]/1000)))
    filename = filename_base + x_coord +'_'+ y_coord + '.nc'

    terminal_inversion_smooth(m, C, p_b, p_c, erB, erC, n, adj, square_size, tapering, centre_include, \
                              data[i], trans_funcs, \
                                  filepath_itslive, filepath_rema, filepath_bedmach, interp_grid_spacing, \
                                  CutB, CutC, wavcutC, wavcutB, filename) 
    
#else:

#file_prefix = 'PIG_drot_REMA_29_16_'
#filename = file_prefix + str(rank) + '.nc'
#print(filename)

