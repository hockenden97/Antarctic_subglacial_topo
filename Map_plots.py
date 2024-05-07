# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:16:11 2023

@author: s1502187
"""

print('Welcome to map_plots.py')

# Load useful packages 
import numpy as np
from inversion_module_v2 import smooth_data_load
from netCDF4 import Dataset
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import box

# Set filepaths to useful data sources
# Datafiles
filepath_rema = '../Data/GaplessREMA100.nc'
filepath_bedmach = "../Data/BedMachineAntarctica_2019-11-05_v01.nc"
filepath_bedmach = "../Data/BedMachineAntarctica-v3.nc"
filepath_itslive = "../Data/antarctica_ice_velocity_450m_v2.nc"
filepath_ifpa = "IFPA_bed.nc"
groundinglineM = gpd.read_file("GroundingLine_Antarctica_v2.shp")

bounds_ar = [[-1.82e6, -1.45e6, -0.25e6, 0.075e6], 
            [-1.45e6, -1.0e6, -0.5e6,0e6], 
            [-1.25e6, -0.75e6, -1.0e6, -0.5e6],            
            [-1.0e6, -0.55e6, -0.5e6,0e6],  
            [-0.55e6, 0e6, -0.4e6, 0.1e6], 
            [-0.55e6, 0e6, 0.1e6, 0.6e6], 
            [-0.5e6, 0e6, 0.6e6, 1.1e6], 
            [-0.5e6, 0e6, 1.1e6, 1.6e6], 
            [0e6, 0.5e6, -0.6e6, -0.1e6], 
            [0e6, 0.5e6, -0.1e6, 0.4e6], 
            [0e6, 0.5e6, 0.4e6, 0.9e6], 
            [0e6, 0.5e6, 0.9e6, 1.4e6], 
            [0e6, 0.5e6, 1.4e6, 1.9e6], 
            [0.6e6, 1.0e6, -1.9e6, -1.5e6],
            [0.6e6, 1.0e6, -1.5e6, -1.0e6],
            [0.5e6, 1.0e6, -1.0e6, -0.5e6], 
            [0.5e6, 1.0e6, -0.5e6, 0.0e6], 
            [0.5e6, 1.0e6, 0.0e6, 0.5e6], 
            [0.5e6, 1.0e6, 0.5e6, 1.0e6], 
            [0.5e6, 1.0e6, 1.0e6, 1.5e6], 
            [0.5e6, 1.0e6, 1.5e6, 1.8e6], 
            [1.0e6, 1.45e6, -1.9e6, -1.47e6], 
            [1.0e6, 1.5e6, -1.47e6, -1.0e6], 
            [1.0e6, 1.5e6, -1.0e6, -0.5e6], 
            [1.0e6, 1.5e6, -0.5e6, 0.0e6], 
            [1.0e6, 1.5e6, 0.0e6, 0.5e6], 
            [1.0e6, 1.5e6, 0.5e6, 1.0e6], 
            [1.0e6, 1.5e6, 1.0e6, 1.5e6], 
            [1.45e6, 1.75e6, -1.75e6, -1.47e6], 
            [1.5e6, 1.95e6, -1.47e6, -0.97e6], 
            [1.5e6, 1.9e6, -0.97e6, -0.5e6], 
            [1.9e6, 2.3e6, -0.97e6, -0.5e6], 
            [1.5e6, 2.0e6, -0.5e6, 0.0e6], 
            [2.0e6, 2.4e6, -0.5e6, 0.0e6], 
            [1.5e6, 2.0e6, -0.0e6, 0.5e6], 
            [1.5e6, 1.95e6, 0.9e6, 1.4e6]]

clims = [[-1300,800], 
        [-2000,1000], 
        [-1500,1500], 
        [-2000,1500], 
        [-1000,2500], 
        [-1500,1500], 
        [-2000,1500], 
        [-1500,1300],
        [-600,2900],
        [-1000,1000],
        [-1000,1000],
        [-700,700],
        [-300,900],
        [-1200,1600],
        [-900,1100],
        [-700,700],
        [-400,1300],
        [-200,2000],
        [-500,1200],
        [0,1600],
        [-500,2500],
        [-1700,700],
        [-1200,800],
        [-1200,800],
        [-1200,1500],
        [-1000,2500],
        [-1000,1500],
        [-200,1400],
        [-1500,600],
        [-1100,600],
        [-1500,400],
        [-1600,500],
        [-600,1200],
        [-1400,1000],
        [-500,1400],
        [-500,1600]]


for i in range(len(bounds_ar)):
#for i in (0,1):
    # Extend bounds slightly for importing data
    bounds_import = [bounds_ar[i][0] - 200, bounds_ar[i][1] + 200, bounds_ar[i][2] - 200, bounds_ar[i][3] + 200]
    # Import the ITSLIVE, REMA and Bedmachine data
    X_its, Y_its, VX, VY, X_rema, Y_rema, SURF, X_bedmach, Y_bedmach, thick, bedmach, errbed, source = \
    smooth_data_load(bounds_import, filepath_itslive, filepath_rema, filepath_bedmach)
    # To load the IFPA data
    fh = Dataset(filepath_ifpa, 'r', format = 'NETCDF4');
    X = fh.variables['x'][:]
    Y = fh.variables['y'][:]
    xl = next(x for x, val in enumerate(X) if val >= bounds_import[0]) # X in ascending order
    xh = next(x for x, val in enumerate(X) if val >= bounds_import[1])
    yl = next(x for x, val in enumerate(Y) if val >= bounds_import[2])  # Y in descending order
    yh = next(x for x, val in enumerate(Y) if val >= bounds_import[3])  
    bed_ifpa = fh.variables['bed'][yl:yh, xl:xh]
    X_ifpa = fh.variables['x'][xl:xh]
    Y_ifpa = fh.variables['y'][yl:yh]
    X_ifpa, Y_ifpa = np.meshgrid(np.array(X_ifpa), np.array(Y_ifpa))
    fh.close()
    
    ax = [[],[],[],[]]
    fig = plt.figure(constrained_layout=True, figsize = (20,8))
    widths = [6,0.01,6,4]
    gs = fig.add_gridspec(1,4, width_ratios = widths)
    ax[0] = fig.add_subplot(gs[0, 0])
    ax[1] = fig.add_subplot(gs[0, 1])
    ax[2] = fig.add_subplot(gs[0, 2])
    ax[3] = fig.add_subplot(gs[0, 3])
    
    d = 5
    d2 = 1
    im = [[],[]]
    im[0] = ax[0].pcolor(X_bedmach[::d2,::d2], Y_bedmach[::d2,::d2], bedmach[::d2,::d2])
    im[1] = ax[2].pcolor(X_ifpa[::d,::d], Y_ifpa[::d,::d], bed_ifpa[::d,::d])
    #cbar2  = plt.colorbar(im[0], cax = ax[1], shrink = 0.5)
    titles = ('\nBedmachine Antarctica', '\nIce Flow\nPerturbation Analysis')
    #clims = (np.min((np.min(bedmach), np.nanmin(bed_ifpa))), np.max((np.max(bedmach), np.nanmax(bed_ifpa))))
    
    for j in (0,2, 3):
        groundinglineM.plot(ax = ax[j], facecolor = 'None', edgecolor='k')
        ax[j].get_xaxis().set_ticks([]); ax[j].get_yaxis().set_ticks([]);   
        
    
    for j in range(2):
        ax[2*j].set_xlim(bounds_ar[i][0], bounds_ar[i][1])
        ax[2*j].set_ylim(bounds_ar[i][2], bounds_ar[i][3])
        ax[2*j].set_title(titles[j], fontsize = 30)
        im[j].set_clim(clims[i])
   
    cbar = plt.colorbar(im[0], ax = ax[0], shrink = 0.5)
    cbar.set_label('Elevation (m)', fontsize = 20)
    cbar.ax.tick_params(labelsize=16)
    ax[1].axis('off');
    ax[3].axis('off');
    
    x_width  = (bounds_ar[i][1] - bounds_ar[i][0])/1000
    y_width =  (bounds_ar[i][3] - bounds_ar[i][2])/1000
    cbar.ax.annotate('Area dimensions\n{:.0f} km by {:.0f} km\n'.format(x_width, y_width), \
                   xy=(2, 1), xytext=(0, 1.02), xycoords = 'axes fraction', fontsize = 16);
    
    for kk in range(len(bounds_ar)):
        geom = box(bounds_ar[kk][0], bounds_ar[kk][2], bounds_ar[kk][1], bounds_ar[kk][3])
        ax[3].plot(*geom.exterior.xy, color = '#bebebe')
        ax[3].annotate(kk + 1, xy=(2, 1), xytext=(np.mean((bounds_ar[kk][0], bounds_ar[kk][1])), \
                                              np.mean((bounds_ar[kk][2], bounds_ar[kk][3]))), \
                    fontsize = 12, ha = 'center', va = 'center', color = '#bebebe');
    geom = box(bounds_ar[i][0], bounds_ar[i][2], bounds_ar[i][1], bounds_ar[i][3])
    ax[3].plot(*geom.exterior.xy, color = 'r')
    ax[3].annotate(i + 1, xy=(2, 1), xytext=(np.mean((bounds_ar[i][0], bounds_ar[i][1])), \
                                              np.mean((bounds_ar[i][2], bounds_ar[i][3]))), \
                    fontsize = 12, ha = 'center', va = 'center', color = 'r');
    
    #plt.title('Canyon Basin', fontsize = 20)
    
    fig.savefig('../Nature_figures/Extra_topo/Extra_topo_{}.png'.format(i+1), dpi = 200, bbox_inches = 'tight')
    plt.close(fig)
    print(i, 'out of', len(bounds_ar), 'Plotted', end = '\r')