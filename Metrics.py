"""
This code was created by Helen Ockenden (hockenden97) on the 28th of November 2023.
It was last updated on the 28th of November 2023.

This code takes the IFPA tiles and the Bedmachine topography and calculates a selection of statistical measures
- Mean elevation
- Standard deviation
- Standard deviation with slope removed
- Skewness
- Kurtosis
- RMS slope
- RMS curvature
- Low and high frequency standard deviation
- Low and high frequency RMS slope
- Low and high frequency RMS curvature
- Fractal dimension (wavelengths >5km)
- Fractal dimension (wavelengths > ice thickness)
- Wavelength of maximum spectral power
- Angle of maximum spectral power
- The number of 20m hills in 5km
- The number of 50m hills in 5km
- The number of 100m hills in 5km
- The number of 250m hills in 5km

These metrics are saved to .nc
"""

# Load in useful packages 
import numpy as np
import glob
from netCDF4 import Dataset
from inversion_module_v2 import smooth_data_load
from skimage.filters import difference_of_gaussians
from scipy.fft import fft2
from sklearn.linear_model import LinearRegression
from scipy.stats import skew, kurtosis 
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def slope_removal(grid, return_angle = False):
    y_grid = np.arange(0,grid.shape[0], 1)
    y_grid = y_grid - np.mean(y_grid)
    x_grid = np.arange(0,grid.shape[1], 1)
    x_grid = x_grid - np.mean(x_grid)
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    xs = np.ndarray.flatten(x_grid[::,::])
    ys = np.ndarray.flatten(y_grid[::,::])
    zs = np.ndarray.flatten(grid[::, ::])
    # do fit
    tmp_A = []
    tmp_b = []
    for i in range(len(xs)):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_b.append(zs[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    # Manual solution
    fit = (A.T * A).I * A.T * b
    #alpha_s = np.array((fit[0]**2 + fit[1] **2)/ (np.sqrt((fit[0]**2 + fit[1] **2))))
    angle = (np.arctan2(fit[1],fit[0]) *180/np.pi)
    #Calculate the slope array
    Z2 = x_grid * np.array(fit[0]) + y_grid * np.array(fit[1]) + np.array(fit[2])
    if return_angle == True:
        return angle
    else:
        return Z2

# Create a function to calculate the desired hypsometric parameters
def bandpass_analysis(grid):
    # Remove the mean slope of the bed
    bed_deslope = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) 
    # Set filtering parameters
    filtering_params = (4,80) # High frequency, low frequency
    # Create variables to fill
    rms_slope = np.ndarray(len(filtering_params))
    rms_curv = np.ndarray(len(filtering_params))
    std = np.ndarray(len(filtering_params))
    for i in range(len(filtering_params)):
        # Bandpass filter the images
        filtered_images = difference_of_gaussians(bed_deslope, filtering_params[i])
        # Calculate slope
        slope = np.gradient(filtered_images)
        rms_slope[i] = np.sqrt(np.mean(slope[0]**2 + slope[1] **2))
        # Calculate curvature
        curvature = np.gradient(slope)
        rms_curv[i] = np.sqrt(np.mean(curvature[0]**2 + curvature[1] **2))
        # Calculate standard deviation
        std[i] = np.std(filtered_images[i])
    return rms_slope[0], rms_curv[0], std[0], \
           rms_slope[1], rms_curv[1], std[1]

# Create a function to do Hann tapering 
def tapering_hann(a, tapering):
    x_taper = int(a.shape[1])
    y_taper = int(a.shape[0])
    x_sin = (np.sin(np.linspace(-np.pi/2,3*np.pi/2,x_taper)) +1 ) /2
    y_sin = (np.sin(np.linspace(-np.pi/2,3*np.pi/2,y_taper)) +1 ) /2
    x_sin, y_sin = np.meshgrid(x_sin, y_sin)
    z = x_sin * y_sin
    return z

# Create a function to do a fourier transform
def fourier_transform(grid, x_grid, y_grid):
    spacing = np.abs(x_grid[0,0] - x_grid[0,1])
    grid_for_fft = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) *tapering_hann(grid, 0.1)
    grid_for_fft = grid_for_fft - np.mean(grid_for_fft)
    grid_fft2 = fft2(grid_for_fft)
    ar1 = np.fft.fftfreq(grid.shape[1], spacing)
    ar2 = np.fft.fftfreq(grid.shape[0], spacing)
    k,l = np.meshgrid(ar1,ar2)
    j = np.sqrt(k**2 + l**2)
    return grid_fft2, j 

# Create a function to calculate the power spectral density 
def PSD_smoothed(grid, x_grid, y_grid, filtering):
    spacing = np.abs(x_grid[0,0] - x_grid[0,1])
    grid_fft2, j = fourier_transform(grid, x_grid, y_grid)
    j[0,0] = np.nan
    Nx = grid.shape[1]
    Ny = grid.shape[0]
    W_sq_sum = np.sum(tapering_hann(grid, 0.1) ** 2)
    #P_dft_w = np.abs(grid_fft2) **2 / (Nx * Ny * W_sq_sum)
    #P_psd = P_dft_w * Nx * Ny * spacing * spacing
    # Bandpass filter
    mask = (1/j) > filtering # wavelength in m
    grid_fft2[~mask] = 0
    P_dft_w_smoothed = np.abs(grid_fft2) **2 / (Nx * Ny * W_sq_sum)
    P_psd_smooth = P_dft_w_smoothed * Nx * Ny * spacing * spacing
    # Calculate linear slope 
    j_model  = np.log10(j[mask]).reshape(-1,1)
    p_model = np.log10(P_psd_smooth[mask])
    model = LinearRegression().fit(j_model, p_model)
    b1_roll = model.coef_ 
    return b1_roll

# Create a second function to do power spectral density and return the spectrum
def power_spectral_density(grid, x_grid, y_grid):
    spacing = np.abs(x_grid[0,0] - x_grid[0,1])
    grid_for_fft = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) *tapering_hann(grid, 0.1)
    grid_for_fft = grid_for_fft - np.mean(grid_for_fft)
    grid_fft2 = fft2(grid_for_fft)
    Nx = grid.shape[1]
    Ny = grid.shape[0]
    P_dft = np.abs(grid_fft2) **2 / (Nx ** 2 * Ny **2)
    P_psd = P_dft * Nx * Ny * spacing * spacing
    ar1 = np.fft.fftfreq(grid.shape[1], spacing)
    ar2 = np.fft.fftfreq(grid.shape[0], spacing)
    k,l = np.meshgrid(ar1,ar2)
    j = np.sqrt(k**2 + l**2)
    theta = np.arctan2(l,k)
    return P_psd, j, theta

# Fit a selection of regression lines for various frequencies
def frequency_regression(P_psd, j, theta, plotting = False):
    j_flat = np.log10(np.ndarray.flatten(j)[1:])
    no_j_values = 72
    j_regression = (np.linspace(np.nanmin(j_flat), np.nanmax(j_flat), no_j_values))
    P_psd_flat = np.log10(np.ndarray.flatten(P_psd))[1:]
    b0 = np.ndarray(len(j_regression))
    b1 = np.ndarray(len(j_regression))
    r_sq = np.ndarray(len(j_regression))
    for i in range(len(j_regression)):
        mask = j_flat <= j_regression[i]
        j_flat_reg = j_flat[mask].reshape((-1, 1))
        P_psd_flat_reg = P_psd_flat[mask]
        model = LinearRegression().fit(j_flat_reg, P_psd_flat_reg)
        b0[i], b1[i], r_sq[i] = model.intercept_ , model.coef_ , model.score(j_flat_reg, P_psd_flat_reg)
    P_psd_regression = j_regression * b1 + b0   
    grad_b1  = np.gradient(np.abs(b1))
    mask = j_regression > -4.0
    for i in range(len(j_regression[mask])):
        if grad_b1[mask][i] == np.max(grad_b1[mask]):
            roll_off = j_regression[mask][i]
    # Creating a line of best fit above the roll-off wavelength    
    mask = j_regression <= roll_off - 0.2    
    j_val_reg = j_regression[mask].reshape(-1,1)
    p_val_reg = P_psd_regression[mask]
    model = LinearRegression().fit(j_val_reg, p_val_reg)
    b0_roll, b1_roll = model.intercept_ , model.coef_ 
  #  p_psd_roll = j_regression * b1_roll + b0_roll
      # Creating a background spectrum (based on the line of best fit)
    j_v2 = j.copy()
    j_v2[0,0] = np.nan
    p_sim = 10 ** (np.log10(j_v2) * b1_roll + b0_roll)
    Corrected_power = P_psd/p_sim
    mask = np.log10(j_v2) > roll_off
    Corrected_power[mask] = 0
    Corrected_power[0,0] = 0
    # What are the values of theta and j at the maximum corrected power?
    mask = Corrected_power == np.max(Corrected_power)
    j_max_power, theta_max_power = j[mask], theta[mask]
 #   if plotting == True:
 #       return P_psd_regression, p_psd_roll, j_regression, 
    return 1/j_max_power[0], theta_max_power[0] * 180/np.pi  

# Create a function to calculate mean, std, slope and curvature 
def hypso2_analysis(bed_topo_dem):
    mean = np.mean(bed_topo_dem)
    std = np.std(bed_topo_dem)
    grid = bed_topo_dem
    bed_deslope = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) 
    std_deslope = np.std(bed_deslope)
    dem = bed_topo_dem
    x, y = np.gradient(dem)
    slope = np.arctan(np.sqrt(x**2 + y**2))
    rms_slope = np.sqrt(np.mean(slope **2))
    x2, y2 = np.gradient(slope)
    curvature = np.arctan(np.sqrt(x2**2 + y2 **2))
    rms_curvature = np.sqrt(np.mean(curvature **2))
    return mean, std, std_deslope, rms_slope, rms_curvature

# Create a function to calculate relief, skewness and kurtosis
def hypso_analysis(bed_topo_dem):
    """ A function to calculate various hypsometric parameters for the provided DEM """
    """ Relief = Max_elevation - Min_elevation                                      """
    """ Skewness, from the scipy.stats.skew function                                """
    """ Kurtosis, from the scipy.stats.kurtosis function                            """
    # Flatten the 2D grid to get a single distribution
    bed_topo = np.ndarray.flatten(bed_topo_dem)
    relief = np.max(bed_topo) - np.min(bed_topo)
    bed_skewness = skew(np.ndarray.flatten(bed_topo))
    bed_kurt = kurtosis(np.ndarray.flatten(bed_topo))
    return relief, bed_skewness, bed_kurt

# Create a function to count the number of bumps of a certain threshold 
# Neighbourhood size is baked in as 5km here but could be changed
def bumpiness(bed_ifpa, bedmach, threshold_a):
    grid = bed_ifpa
    bed_deslope = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) 
    grid = bedmach
    bedmach_deslope = ((grid - slope_removal(grid))-np.mean(grid-slope_removal(grid))) 
    data = bed_deslope, bedmach_deslope
    neighborhood_size = 50, 10
    threshold = threshold_a, threshold_a 
    ifpa_count = []
    bedmach_count = []
    for j in range(len(data)):
        data_max = filters.maximum_filter(data[j], neighborhood_size[j])
        maxima = (data[j] == data_max)
        data_min = filters.minimum_filter(data[j], neighborhood_size[j])
        diff = ((data_max - data_min) > threshold[j])
        maxima[diff == 0] = 0

        labeled, num_objects = ndimage.label(maxima)
        slices = ndimage.find_objects(labeled)
        x, y = [], []
        for dy,dx in slices:
            x_center = (dx.start + dx.stop - 1)/2
            y_center = (dy.start + dy.stop - 1)/2
            # To avoid double counting between squares we remove maxima on the 0 edges 
            if x_center == 0:
                None
            else: 
                if y_center == 0:
                    None
                else:
                    x.append(x_center)
                    y.append(y_center)
        if j == 0:
                ifpa_count = len(x)
        else:
            bedmach_count = len(x)
    return ifpa_count, bedmach_count

           
# Load in the list of REMA grid filenames 
filepaths = glob.glob('../Inversion_code/*A_grids/*A*.nc')
# Update terminal on progress
print('There are', len(filepaths), 'tiles')

# Set the filepath for the Bedmachine data (should this be v3?)
filepath_rema = '../Data/GaplessREMA100.nc'
filepath_bedmach = "../Data/BedMachineAntarctica_2019-11-05_v01.nc"
filepath_itslive = "../Data/antarctica_ice_velocity_450m_v2.nc"
# Update terminal on progress
print('Setting filepaths to data')

# Create a load of empty variables to save parameters to
# Coordinates
x_ifpa                  = np.ndarray(len(filepaths))
y_ifpa                  = np.ndarray(len(filepaths)) 
# Elevation/Hypsometry
ifpa_mean               = np.ndarray(len(filepaths))
ifpa_std                = np.ndarray(len(filepaths))
ifpa_std_deslope        = np.ndarray(len(filepaths))
ifpa_rms_slope          = np.ndarray(len(filepaths))
ifpa_rms_curvature      = np.ndarray(len(filepaths))
ifpa_relief             = np.ndarray(len(filepaths))
ifpa_elevskew           = np.ndarray(len(filepaths))
ifpa_kurtosis           = np.ndarray(len(filepaths))
bedmach_mean            = np.ndarray(len(filepaths))
bedmach_std             = np.ndarray(len(filepaths))
bedmach_std_deslope     = np.ndarray(len(filepaths))
bedmach_rms_slope       = np.ndarray(len(filepaths))
bedmach_rms_curvature   = np.ndarray(len(filepaths))
bedmach_relief          = np.ndarray(len(filepaths))
bedmach_elevskew        = np.ndarray(len(filepaths))
bedmach_kurtosis        = np.ndarray(len(filepaths))
# Bandpass filtered topography
ifpa_rms_slope_h        = np.ndarray(len(filepaths))
ifpa_rms_curv_h         = np.ndarray(len(filepaths))
ifpa_std_h              = np.ndarray(len(filepaths))
ifpa_rms_slope_l        = np.ndarray(len(filepaths))
ifpa_rms_curv_l         = np.ndarray(len(filepaths))
ifpa_std_l              = np.ndarray(len(filepaths))
bedmach_rms_slope_h     = np.ndarray(len(filepaths))
bedmach_rms_curv_h      = np.ndarray(len(filepaths))
bedmach_std_h           = np.ndarray(len(filepaths))
bedmach_rms_slope_l     = np.ndarray(len(filepaths))
bedmach_rms_curv_l      = np.ndarray(len(filepaths))
bedmach_std_l           = np.ndarray(len(filepaths))
# Spectral properties
ifpa_b1_5km             = np.ndarray(len(filepaths))
ifpa_b1_thickness       = np.ndarray(len(filepaths))
ifpa_wav_max_power      = np.ndarray(len(filepaths))
ifpa_theta_max_power    = np.ndarray(len(filepaths))
bedmach_b1_5km          = np.ndarray(len(filepaths))
bedmach_b1_thickness    = np.ndarray(len(filepaths)) 
bedmach_wav_max_power   = np.ndarray(len(filepaths))
bedmach_theta_max_power = np.ndarray(len(filepaths))
# Bumpiness 
ifpa_count_max_20       = np.ndarray(len(filepaths))
ifpa_count_max_50       = np.ndarray(len(filepaths))
ifpa_count_max_100      = np.ndarray(len(filepaths))
ifpa_count_max_250      = np.ndarray(len(filepaths))
bedmach_count_max_20    = np.ndarray(len(filepaths))
bedmach_count_max_50    = np.ndarray(len(filepaths))
bedmach_count_max_100   = np.ndarray(len(filepaths))
bedmach_count_max_250   = np.ndarray(len(filepaths))

# For each grid in the range desired 
#for i in range(len(filepaths)):
for i in range(5):
    # Load in the IFPA data 
    output = Dataset(filepaths[i],'r',format='netCDF4')
    X = output.variables['x'][:]
    Y = output.variables['y'][:]
    x_i, y_i = np.meshgrid(X, Y)
    x_ifpa[i] = np.mean(x_i)
    y_ifpa[i] = np.mean(y_i)
    bed_i = output.variables['bed'][:,:]
    output.close()
    bounds = [np.min(x_i), np.max(x_i), np.min(y_i), np.max(y_i)]
    # Load in the bedmachine data
    #X_its, Y_its, VX, VY, X_rema, Y_rema, SURF, X_bedmach, Y_bedmach, thick, bedmach, errbed, source = \
    #    smooth_data_load(bounds, filepath_itslive, filepath_rema, filepath_bedmach)
    _, _, _, _, _, _, _, X_bedmach, Y_bedmach, thick, bedmach, _, _ = \
        smooth_data_load(bounds, filepath_itslive, filepath_rema, filepath_bedmach)
    # Elevation/hypsometry parameters
    ifpa_mean[i], ifpa_std[i], ifpa_std_deslope[i], ifpa_rms_slope[i], \
        ifpa_rms_curvature[i] = hypso2_analysis(bed_i)
    bedmach_mean[i], bedmach_std[i], bedmach_std_deslope[i], bedmach_rms_slope[i], \
        bedmach_rms_curvature[i] = hypso2_analysis(bedmach)
    ifpa_relief[i], ifpa_elevskew[i], ifpa_kurtosis[i] = hypso_analysis(bed_i)
    bedmach_relief[i], bedmach_elevskew[i], bedmach_kurtosis[i] = hypso_analysis(bedmach)
    # Calculate the low and high frequency metrics for each bed topography 
    ifpa_rms_slope_h[i], ifpa_rms_curv_h[i], ifpa_std_h[i], \
    ifpa_rms_slope_l[i], ifpa_rms_curv_l[i], ifpa_std_l[i] = \
           bandpass_analysis(bed_i)
    bedmach_rms_slope_h[i], bedmach_rms_curv_h[i], bedmach_std_h[i], \
    bedmach_rms_slope_l[i], bedmach_rms_curv_l[i], bedmach_std_l[i] = \
           bandpass_analysis(bedmach)
    # Calculate the fractal dimension metrics for each bed topography 
    ifpa_b1_5km[i] = PSD_smoothed(bed_i, x_i, y_i, 5000)
    ifpa_b1_thickness[i] = PSD_smoothed(bed_i, x_i, y_i, np.nanmean(thick))
    bedmach_b1_5km[i] = PSD_smoothed(bedmach, X_bedmach, Y_bedmach, 5000)
    bedmach_b1_thickness[i] = PSD_smoothed(bedmach, X_bedmach, Y_bedmach, np.nanmean(thick))
    P_psd_bedmach, j_bedmach, theta_bedmach = power_spectral_density(bedmach, X_bedmach, Y_bedmach)
    P_psd_ifpa, j_ifpa, theta_ifpa = power_spectral_density(bed_i, x_i, y_i)
    bedmach_wav_max_power[i], bedmach_theta_max_power[i] = \
            frequency_regression(P_psd_bedmach, j_bedmach, theta_bedmach)
    ifpa_wav_max_power[i], ifpa_theta_max_power[i] = \
            frequency_regression(P_psd_ifpa, j_ifpa, theta_ifpa)
    # Calculate the bumpiness metrics for each bed topography 
    ifpa_count_max_20[i], bedmach_count_max_20[i] = bumpiness(bed_i, bedmach, 20)
    ifpa_count_max_50[i], bedmach_count_max_50[i] = bumpiness(bed_i, bedmach, 50)
    ifpa_count_max_100[i], bedmach_count_max_100[i] = bumpiness(bed_i, bedmach, 100)
    ifpa_count_max_250[i], bedmach_count_max_250[i] = bumpiness(bed_i, bedmach, 250)
    # Update terminal on progress
    print(i+1, 'out of', len(filepaths), 'processed', end = '\r')

    
# Create a function to save this data to file
def NC_data_1D(filename, file):
    ncfile = Dataset(filename, mode = 'w', format = 'NETCDF4')
    # Make data dimensions
    ncfile.createDimension('x', len(file))
    #Create the axes/variable lengths
    data = ncfile.createVariable('data', np.float64, ('x'))
    #Write in the data
    data[:] = file[:]
    ncfile.close()
    
# Save the metrics produced to file
NC_data_1D('x_ifpa.nc', x_ifpa)
NC_data_1D('y_ifpa.nc', y_ifpa)
# Elevation/Hypsometry
NC_data_1D('ifpa_mean.nc', ifpa_mean)
NC_data_1D('ifpa_std.nc', ifpa_std)
NC_data_1D('ifpa_std_deslope.nc', ifpa_std_deslope)
NC_data_1D('ifpa_rms_slope.nc', ifpa_rms_slope)
NC_data_1D('ifpa_rms_curvature.nc', ifpa_rms_curvature)
NC_data_1D('ifpa_relief.nc', ifpa_relief)
NC_data_1D('ifpa_elevskew.nc', ifpa_elevskew)
NC_data_1D('ifpa_kurtosis.nc', ifpa_kurtosis)
NC_data_1D('bedmach_mean.nc', bedmach_mean)
NC_data_1D('bedmach_std.nc', bedmach_std)
NC_data_1D('bedmach_std_deslope.nc', bedmach_std_deslope)
NC_data_1D('bedmach_rms_slope.nc', bedmach_rms_slope)
NC_data_1D('bedmach_rms_curvature.nc', bedmach_rms_curvature)
NC_data_1D('bedmach_relief.nc', bedmach_relief)
NC_data_1D('bedmach_elevskew.nc', bedmach_elevskew)
NC_data_1D('bedmach_kurtosis.nc', bedmach_kurtosis)
# Bandpass filtered topography
NC_data_1D('i_rms_slope_h.nc', ifpa_rms_slope_h)
NC_data_1D('i_rms_curv_h.nc', ifpa_rms_curv_h) 
NC_data_1D('i_std_h.nc', ifpa_std_h)
NC_data_1D('i_rms_slope_l.nc', ifpa_rms_slope_l) 
NC_data_1D('i_rms_curv_l.nc', ifpa_rms_curv_l) 
NC_data_1D('i_std_l.nc', ifpa_std_l)
NC_data_1D('b_rms_slope_h.nc', bedmach_rms_slope_h) 
NC_data_1D('b_rms_curv_h.nc', bedmach_rms_curv_h) 
NC_data_1D('b_std_h.nc', bedmach_std_h) 
NC_data_1D('b_rms_slope_l.nc', bedmach_rms_slope_l) 
NC_data_1D('b_rms_curv_l.nc', bedmach_rms_curv_l) 
NC_data_1D('b_std_l.nc', bedmach_std_l)
# Spectral characteristics
NC_data_1D('ifpa_b1_5km.nc', ifpa_b1_5km)
NC_data_1D('ifpa_b1_thickness.nc', ifpa_b1_thickness)
NC_data_1D('ifpa_wav_max_power.nc', ifpa_wav_max_power)
NC_data_1D('ifpa_theta_max_power.nc', ifpa_theta_max_power)
NC_data_1D('bedmach_wav_max_power.nc', bedmach_wav_max_power)
NC_data_1D('bedmach_theta_max_power.nc', bedmach_theta_max_power)
NC_data_1D('bedmach_b1_5km.nc', bedmach_b1_5km)
NC_data_1D('bedmach_b1_thickness.nc', bedmach_b1_thickness)
# Bumpiness 
NC_data_1D('ifpa_count_max_20.nc', ifpa_count_max_20)       
NC_data_1D('ifpa_count_max_50.nc', ifpa_count_max_50)       
NC_data_1D('ifpa_count_max_100.nc', ifpa_count_max_100)       
NC_data_1D('ifpa_count_max_250.nc', ifpa_count_max_250)       
NC_data_1D('bedmach_count_max_20.nc', bedmach_count_max_20)       
NC_data_1D('bedmach_count_max_50.nc', bedmach_count_max_50)       
NC_data_1D('bedmach_count_max_100.nc', bedmach_count_max_100)       
NC_data_1D('bedmach_count_max_250.nc', bedmach_count_max_250)       
# Update terminal on progress
print('\n Data written to file')





