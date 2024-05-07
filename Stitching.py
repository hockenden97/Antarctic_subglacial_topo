from netCDF4 import Dataset
import glob
import geopandas as gpd
from file_processing import bed_coord_NC
import numpy as np
import itertools
from inversion_module_v2 import smooth_data_load
import matplotlib.pyplot as plt

#### Files that may be useful for plotting

groundinglineM = gpd.read_file("GroundingLine_Antarctica_v2.shp")
filepath_rema = '../Data/GaplessREMA100.nc'
filepath_rema_its = '../Data/GaplessREMA100_its2.nc'
filepath_bedmach = "../Data/BedMachineAntarctica-v3.nc"
filepath_itslive_old = "../Data/antarctica_ice_velocity_450m_v2.nc"
filepath_itslive = '../Data/antarctic_ice_vel_phase.nc'

# =============================================================================
# output = Dataset('X_ifpa.nc','r',format='netCDF4')
# #x_ifpa = output.variables['data'][:]
# output.close()
# output = Dataset('Y_ifpa.nc','r',format='netCDF4')
# #y_ifpa = output.variables['data'][:]
# output.close()
# 
# =============================================================================

filepaths = glob.glob('../Dartmouth/C_50/Grids_*.nc')
print('There are', len(filepaths), 'tiles')

# Create a new list of the coordinate points
#x = filepaths[10].split('_')
#x, x[3], x[4]
#x[3].split('m'), '1000'.split('m')
filepath_coords = ([])
for i in range(len(filepaths)):
    breakup = filepaths[i].split('_')
    x_a = breakup[3].split('m')
    y_a = breakup[4].split('m')
    if len(x_a) == 2:
        x = -int(x_a[1]) * 1e3
    else:
        x = int(x_a[0]) * 1e3
    if len(y_a) == 2:
        y = -int(y_a[1].split('.')[0]) * 1e3
    else: 
        y = int(y_a[0].split('.')[0]) * 1e3
    filepath_coords.append([x,y])

np.savetxt('../Dartmouth/C_true/coords_which_ran.csv', filepath_coords, delimiter = ',')

        
fh = Dataset(filepath_itslive, 'r', format='NETCDF4');
X = fh.variables['x'][:]
Y = fh.variables['y'][:]
fh.close()

# Create the central coordinates... 
a = ([0])
for i in range(int(len(X/50000))):
    new_coord = 50000*(i+1) 
    if new_coord +100000 < X.max():
        a.append(new_coord)
    if -new_coord -100000 > X.min():
        a.insert(0,-new_coord)
b = ([0])
for i in range(int(len(Y/50000))):
    new_coord = 50000*(i+1) 
    if new_coord + 100000 < Y.max():
        b.append(new_coord)
    if -new_coord -100000 > Y.min():
        b.insert(0,-new_coord)
        
bounds = [-2.15e6, 2.45e6, -2.0e6, 2.2e6]    
bounds = [-2.0e6, 2.45e6, -2e6, 2.2e6]   
#bounds = [2.1e6,2.4e6, -0.6e6,-0.3e6]
x_cen, y_cen = np.meshgrid(a,b)
x_cen_coords = ([])
y_cen_coords = ([])
for i in range(len(a)):
    if x_cen[0,i] > bounds[0]:
        if x_cen[0,i] < bounds[1] + 1:
            #print(i, x_cen[0,i])
            x_cen_coords.append(i)
for i in range(len(b)):
    if y_cen[i,0] > bounds[2]:
        if y_cen[i,0] < bounds[3] + 1:
            #print(i, y_cen[i,0])
            y_cen_coords.append(i)
x_cen_coords, y_cen_coords = np.meshgrid(x_cen_coords, y_cen_coords)        

print('Creating coordinate array')

list_coords = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
for i,j in itertools.product(range(x_cen_coords.shape[0]), range(x_cen_coords.shape[1])):
    list_coords[i,j] = x_cen[y_cen_coords[i,j], x_cen_coords[i,j]], y_cen[y_cen_coords[i,j], x_cen_coords[i,j]]
        
index_filepaths = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
for i,j in itertools.product(range(x_cen_coords.shape[0]), range(x_cen_coords.shape[1])):
    for k in range(len(filepaths)):
        if list_coords[i,j][0] == int(filepath_coords[k][0]):
            if list_coords[i,j][1] == int(filepath_coords[k][1]):
                #print(i,j,k)
                index_filepaths[i,j] = k
        if list_coords[i,j][0] == 2250000:
            if 2249000 == int(filepath_coords[k][0]):
                if list_coords[i,j][1] == int(filepath_coords[k][1]):
                    #print(i,j,k)
                    index_filepaths[i,j] = k

print('Loading in individual cells')      
x_list = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
y_list = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
bed_list = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
shape = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
bounds2 = np.ndarray(x_cen_coords.shape, dtype = np.ndarray)
for i,j in itertools.product(range(x_cen_coords.shape[0]), range(x_cen_coords.shape[1])):
    bounds = [list_coords[i,j][0] - 25005, list_coords[i,j][0] + 24995, \
              list_coords[i,j][1] - 25001, list_coords[i,j][1] + 24999]
    bounds2[i,j] = bounds
    xl = next(x for x, val in enumerate(X) if val >= bounds[0]) # X in ascending order
    xh = next(x for x, val in enumerate(X) if val >= bounds[1])
    yl = next(x for x, val in enumerate(Y) if val <= bounds[2])  # Y in descending order
    yh = next(x for x, val in enumerate(Y) if val <= bounds[3])  
    X_its = X[xl:xh]
    Y_its = np.flip(Y[yh:yl],0)   
    x_all, y_all = np.meshgrid(X_its, Y_its)
    #print(i,j, X_its.shape, Y_its.shape)
    if index_filepaths[i,j] == None:
        #x_fake = np.linspace(list_coords[i,j][0] - 24950, list_coords[i,j][0] + 24950, 500)
        #y_fake = np.linspace(list_coords[i,j][1] - 24950, list_coords[i,j][1] + 24950, 500)
        x_fake = X_its
        y_fake = Y_its
        x_list[i,j], y_list[i,j] = np.meshgrid(x_fake, y_fake)
        bed_list[i,j] = np.zeros(x_list[i,j].shape) * np.nan
       
    else:
        output = Dataset(filepaths[index_filepaths[i,j]],'r',format='netCDF4')
        X_ifpa = output.variables['x'][:]
        Y_ifpa = output.variables['y'][:]
        x_list[i,j], y_list[i,j] = np.meshgrid(X_ifpa, Y_ifpa)
        bed = output.variables['bed'][:,:]
        if np.sum(np.isnan(bed)) != 0:
            print(i,j,np.sum(np.isnan(bed)), bed.shape[0]*bed.shape[1])
        if bed.shape != x_all.shape:
            #print(i,j, bed.shape, X_its.shape, Y_its.shape)
            #print(bounds)
            #print(np.min(X_ifpa), np.max(X_ifpa), np.min(Y_ifpa), np.max(Y_ifpa))
            #print(np.min(X_its), np.max(X_its), np.min(Y_its), np.max(Y_its)) 
            x_list[i,j], y_list[i,j] = np.meshgrid(X_ifpa[:len(X_ifpa)-1], Y_ifpa)
            bed = output.variables['bed'][:,:len(X_ifpa)-1]
            #print(bed.shape, X_its.shape, Y_its.shape)
        bed_list[i,j] = bed
        output.close()   
    shape[i,j] = x_list[i,j].shape
    print(i*x_cen_coords.shape[0] + j, 'out of', x_cen_coords.shape[0] * x_cen_coords.shape[1], 'processed', end='\r')


# =============================================================================
# for i,j in itertools.product(range(x_cen_coords.shape[0]), range(x_cen_coords.shape[1])):
#     if bed_list[i,j].shape == (500,500):
#         ijk = 0
#     else:
#         print(i,j, bed_list[i,j].shape)
#         print(index_filepaths[i,j], list_coords[i,j])
# =============================================================================

print('Merging data files')
x_temp = np.ndarray(x_cen_coords.shape[0], dtype = np.ndarray)
y_temp = np.ndarray(x_cen_coords.shape[0], dtype = np.ndarray)
bed_temp = np.ndarray(x_cen_coords.shape[0], dtype = np.ndarray)
for i in range(x_cen_coords.shape[0]):
    x_temp[i], y_temp[i], bed_temp[i] = np.hstack(x_list[i,:]), np.hstack(y_list[i,:]), np.hstack(bed_list[i,:])
x_fin = np.vstack(x_temp[:])
y_fin = np.vstack(y_temp[:])
bed_fin = np.vstack(bed_temp[:])       


print('Saving to file')        
filename = '../Dartmouth/C_true/IFPA_bed_Cvar.nc'
bed_coord_NC(x_fin, y_fin, bed_fin, filename)

bounds = (np.min(x_fin), np.max(x_fin), np.min(y_fin), np.max(y_fin))
_, _, _, _, _, _, _, X_bedmach, Y_bedmach, _, bedmach, _, _ = \
    smooth_data_load(bounds, filepath_itslive, filepath_rema, filepath_bedmach)

div = 20
div2  = 20


print('Generating image')
fig, ax = plt.subplots(1,2, figsize = (20,8))
im = [[],[]]
im[0] = ax[0].pcolor(x_fin[::div, ::div], y_fin[::div, ::div], bed_fin[::div, ::div])
im[1] = ax[1].pcolor(X_bedmach[::div2, ::div2], Y_bedmach[::div2, ::div2], bedmach[::div2, ::div2])
ax[0].set_title('IFPA')
ax[1].set_title('Bedmachine')

for i in range(2):
    ax[i].set_xlim(bounds[0], bounds[1])
    ax[i].set_ylim(bounds[2], bounds[3])
    plt.colorbar(im[i], ax = ax[i], shrink = 0.7)
    groundinglineM.plot(ax = ax[i], facecolor = 'None', edgecolor='k')       
fig.savefig('../Dartmouth/C_true/Comparison_whole_continent_Cvar.png', dpi = 200, bbox_inches = 'tight')
plt.close(fig)    
        
        
        
        
        
        