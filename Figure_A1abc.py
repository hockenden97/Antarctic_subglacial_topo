# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:16:11 2023

@author: s1502187
"""

print('Welcome to map_plot.py')

# Load useful packages 
import numpy as np
from inversion_module_v2 import smooth_data_load
from netCDF4 import Dataset
import geopandas as gpd
import matplotlib.pyplot as plt
import itertools

# Set filepaths to useful data sources
# Datafiles
filepath_rema = '../Data/GaplessREMA100.nc'
filepath_bedmach = "../Data/BedMachineAntarctica_2019-11-05_v01.nc"
filepath_bedmach = "../Data/BedMachineAntarctica-v3.nc"
filepath_itslive = "../Data/antarctica_ice_velocity_450m_v2.nc"
filepath_ifpa = "IFPA_bed.nc"
groundinglineM = gpd.read_file("GroundingLine_Antarctica_v2.shp")


#bounds = [-1.975e6, 2.425e6, -1.475e6, -0.975e6]

# Extend bounds slightly for importing data
#bounds_import = [bounds[0] - 200, bounds[1] + 200, bounds[2] - 200, bounds[3] + 200]
# To load the IFPA data
fh = Dataset(filepath_ifpa, 'r', format = 'NETCDF4');
#X = fh.variables['x'][:]
#Y = fh.variables['y'][:]
#xl = next(x for x, val in enumerate(X) if val >= bounds_import[0]) # X in ascending order
#xh = next(x for x, val in enumerate(X) if val >= bounds_import[1])
#yl = next(x for x, val in enumerate(Y) if val >= bounds_import[2])  # Y in descending order
#yh = next(x for x, val in enumerate(Y) if val >= bounds_import[3])  
bed_ifpa = fh.variables['bed'][:,:]#[yl:yh, xl:xh]
X_ifpa = fh.variables['x'][:]#[xl:xh]
Y_ifpa = fh.variables['y'][:]#[yl:yh]
X_ifpa, Y_ifpa = np.meshgrid(np.array(X_ifpa), np.array(Y_ifpa))
fh.close()
print('Data loaded')

fh = Dataset(filepath_bedmach, 'r', format='NETCDF4');
X = fh.variables['x'][:]
Y = fh.variables['y'][:]
bed = fh.variables['bed'][:,:]
X_bedmach, Y_bedmach = np.meshgrid(X, Y)
fh.close()

ax = [[],[],[],[],[]]
fig = plt.figure(constrained_layout=True, figsize = (14,5))
widths = [2,2,2,2,2,2]
heights = [4,0.2]
gs = fig.add_gridspec(2,6, width_ratios = widths, height_ratios = heights)
ax[0] = fig.add_subplot(gs[0,0:2])
ax[1] = fig.add_subplot(gs[0,2:4])
ax[2] = fig.add_subplot(gs[0,4:6])
ax[3] = fig.add_subplot(gs[1, 1:3])
ax[4] = fig.add_subplot(gs[1, 3:5])

d = 20
clims = [-2000,2000]
im = ax[0].pcolor(X_ifpa[::d,::d], Y_ifpa[::d,::d], bed_ifpa[::d,::d], zorder = 1)
im.set_clim(clims)
groundinglineM.plot(ax = ax[1], facecolor = '#1d0148', edgecolor='grey', zorder = 0)
d = 25
mask = bed < -2000
bed[mask] = np.nan
im = ax[1].pcolor(X_bedmach[::d,::d], Y_bedmach[::d,::d], bed[::d,::d], zorder = 1)
im.set_clim(clims)
for i in range(3):
    ax[i].set_xlim(-2.6e6, 2.7e6)
    ax[i].set_ylim(-2.3e6, 2.3e6)
ax[4].axis('off')

cb = plt.colorbar(im, cax = ax[3], orientation  = 'horizontal', shrink = 0.5)
cb.set_label('Bed elevation (m)', fontsize = 20, fontweight = 'bold')
cb.ax.tick_params(labelsize=16)
#cb = plt.colorbar(im, cax = ax[4], orientation  = 'horizontal', shrink = 0.5)
#cb.set_label('Bed elevation (m)', fontsize = 20, fontweight = 'bold')
#cb.ax.tick_params(labelsize=16)

mid_letters = ['a)','b)','c)']
for i in range(3):
    groundinglineM.plot(ax = ax[i], facecolor = 'None', edgecolor='grey', zorder = 2)
    ax[i].get_xaxis().set_ticks([]); ax[i].get_yaxis().set_ticks([]);
    ax[i].axis('off')
    ax[i].annotate(mid_letters[i], \
               xy=(2, 1), xytext=(0.09, 0.87), xycoords = 'axes fraction', fontsize = 20);

ax[0].annotate('Ice Flow\nPerturbation\nAnalysis', xy = (2,1), xytext = (0.3, 0.08), ha = 'center', va = 'center', \
                  fontsize = 20, xycoords = 'axes fraction')
ax[1].annotate('Bedmachine\nAntarctica', xy = (2,1), xytext = (0.23, 0.1), ha = 'center', va = 'center', \
                  fontsize = 20, xycoords = 'axes fraction')
fig.tight_layout() 
fig.savefig('../Dartmouth/EXTD1.png', dpi = 200, bbox_inches = 'tight')


# =============================================================================
# fig, ax = plt.subplots(1,1, figsize = (20,20))
# d = 25
# im = ax.pcolor(X_ifpa[::d,::d], Y_ifpa[::d,::d], bed_ifpa[::d,::d])
# groundinglineM.plot(ax = ax, facecolor = 'None', edgecolor='k')
# ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([]); 
# cbar = plt.colorbar(im, ax = ax, shrink = 0.5)
# cbar.set_label('IFPA bed elevation (m)', fontsize = 16)
# im.set_clim(-2000,2000)
# ax.axis('off')
# filename = 'All_Antarctica_map2000.png'
# fig.savefig(filename, dpi = 100, bbox_inches = 'tight')
# plt.close(fig)
# print('Plotting complete and figure saved')
# =============================================================================
