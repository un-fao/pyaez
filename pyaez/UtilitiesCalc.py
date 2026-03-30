"""
PyAEZ version 2.3 (Dec 2024)
Additional calculations used throughout AEZ modules
2020: N. Lakmal Deshapriya
2022/2023: Swun Wunna Htet, K. Boonma
2023 (Dec): Swun Wunna Htet
2024 (Dec): Swun Wunna Htet

Modification:
1. Latitude calculated revised according to GAEZ Fortran routine.
2. New function added: Yield Gap calculation.
3. Removed the python object class feature.
4. Decimal changes for wind speed adjustment on height of measurement based on
   GAEZ routine.
"""

import numpy as np
from scipy.interpolate import interp1d
# try:
#     import gdal
# except:
#     from osgeo import gdal
from osgeo import gdal
np.round_ = np.round
import psutil

def compute_chunk_size_multi(arrays, extra_arrays=None):
    """
    Compute chunk size based on multiple arrays + intermediates.

    Args:
        arrays (list): list of numpy arrays used as input
        extra_arrays (list of dict): intermediate arrays with:
            {
                "shape_factor": tuple relative to (H, W),
                "dtype": dtype
            }
        memory_fraction (float): fraction of available RAM to use
    """
    min_chunk = 128
    max_chunk = 4096

    # Available Memory
    available_mem = psutil.virtual_memory().available
    target_mem = available_mem * 0.1

    # --- Compute bytes per pixel ---
    bytes_per_pixel = 0

    # Input arrays
    for arr in arrays:
        dtype_size = arr.dtype.itemsize

        if arr.ndim == 3:
            depth = arr.shape[2]
            bytes_per_pixel += depth * dtype_size
        else:
            bytes_per_pixel += dtype_size

    # Intermediate arrays
    if extra_arrays:
        for extra in extra_arrays:
            dtype_size = np.dtype(extra["dtype"]).itemsize
            factor = np.prod(extra["shape_factor"])
            bytes_per_pixel += factor * dtype_size

    # Pixels per chunk
    pixels_per_chunk = target_mem // bytes_per_pixel

    # Convert to square chunk
    chunk_dim = int(np.sqrt(pixels_per_chunk))

    # Clamp
    chunk_dim = max(min_chunk, min(chunk_dim, max_chunk))
    print(chunk_dim)

    return chunk_dim

    pixels_per_chunk = target_mem // bytes_per_pixel

    # Convert to square chunk
    chunk_dim = int(np.sqrt(pixels_per_chunk))

    # Clamp
    chunk_dim = max(min_chunk, min(chunk_dim, max_chunk))

    return chunk_dim

def interpMonthlyToDaily( monthly_vector, cycle_begin, cycle_end, no_minus_values=False):
    """Interpolate monthly climate data to daily climate data

    Args:
        monthly_vector (1D NumPy): monthly data that needs interpolating to daily 
        cycle_begin (int): Starting Julian day
        cycle_end (int): Ending Julian day
        no_minus_values (bool, optional): Set minus values to zero. Defaults to False.

    Returns:
        1D NumPy: Daily climate data vector (365 days)
    """        

    doy_middle_of_month = np.arange(0,12)*30 + 15 # Calculate doy of middle of month
    inter_fun = interp1d(doy_middle_of_month, monthly_vector, kind='quadratic', fill_value='extrapolate')
    daily_vector = inter_fun( np.arange(cycle_begin,cycle_end+1) )

    if no_minus_values:
        daily_vector[daily_vector<0] = 0

    return daily_vector

def averageDailyToMonthly(daily_array, leap_year=False,axis=-1):
    """Aggregating daily data into monthly data

    Args:
        daily_vector (1D NumPy Array): daily data array
        leap_year (Boolean): True for leap year, False for non-leap year
    Returns:
        1D NumPy: Monthly data array
    """
    if leap_year:
        month_lengths = np.array([31, 29, 31, 30, 31, 30,
                                  31, 31, 30, 31, 30, 31])
    else:
        month_lengths = np.array([31, 28, 31, 30, 31, 30,
                                  31, 31, 30, 31, 30, 31])

    split_idx = np.cumsum(month_lengths)[:-1]

    daily_array = np.moveaxis(daily_array, axis, -1)

    monthly_sums = np.add.reduceat(daily_array,
                                   np.r_[0, split_idx],
                                   axis=-1)

    monthly_means = monthly_sums / month_lengths

    return monthly_means



    # for i in range(n):
    #     stacked[i] = daily_list[i]
    #
    # # --- allocate output
    # monthly = np.empty(stacked.shape[:-1] + (12,), dtype=stacked.dtype)
    #
    # # --- vectorized monthly aggregation
    # for m in range(12):
    #     start, end = boundaries[m], boundaries[m + 1]
    #     monthly[..., m] = stacked[..., start:end].mean(axis=-1)
    #
    # # --- return as separate arrays
    # return tuple(monthly[i] for i in range(monthly.shape[0]))


    # monthly_vector = np.zeros(12)
    #
    # if leap_year:
    #     monthly_vector[0] = np.sum(daily_vector[:31])/31
    #     monthly_vector[1] = np.sum(daily_vector[31:60])/29
    #     monthly_vector[2] = np.sum(daily_vector[60:91])/31
    #     monthly_vector[3] = np.sum(daily_vector[91:121])/30
    #     monthly_vector[4] = np.sum(daily_vector[121:152])/31
    #     monthly_vector[5] = np.sum(daily_vector[152:182])/30
    #     monthly_vector[6] = np.sum(daily_vector[182:213])/31
    #     monthly_vector[7] = np.sum(daily_vector[213:244])/31
    #     monthly_vector[8] = np.sum(daily_vector[244:274])/30
    #     monthly_vector[9] = np.sum(daily_vector[274:305])/31
    #     monthly_vector[10] = np.sum(daily_vector[305:335])/30
    #     monthly_vector[11] = np.sum(daily_vector[335:])/31
    # else:
    #     monthly_vector[0] = np.sum(daily_vector[:31])/31
    #     monthly_vector[1] = np.sum(daily_vector[31:59])/28
    #     monthly_vector[2] = np.sum(daily_vector[59:90])/31
    #     monthly_vector[3] = np.sum(daily_vector[90:120])/30
    #     monthly_vector[4] = np.sum(daily_vector[120:151])/31
    #     monthly_vector[5] = np.sum(daily_vector[151:181])/30
    #     monthly_vector[6] = np.sum(daily_vector[181:212])/31
    #     monthly_vector[7] = np.sum(daily_vector[212:243])/31
    #     monthly_vector[8] = np.sum(daily_vector[243:273])/30
    #     monthly_vector[9] = np.sum(daily_vector[273:304])/31
    #     monthly_vector[10] = np.sum(daily_vector[304:334])/30
    #     monthly_vector[11] = np.sum(daily_vector[334:])/31
    #
    # return monthly_vector

def generateLatitudeMap(lat_min, lat_max, im_height, im_width):
    """Create latitude map from input geographical extents

    Args:
        lat_min (float): the minimum latitude
        lat_max (float): the maximum latitude
        im_height (float): height of the input raster (pixels,grid cells)
        im_width (float): width of the input raster (pixels,grid cells)

    Returns:
        2D NumPy: interpolated 2D latitude map 
    """        
    lat_step=(lat_max-lat_min)/im_height
    lat_lim = np.linspace(lat_min+lat_step/2, lat_max-lat_step/2, im_height)
    lon_lim = np.linspace(1, 1, im_width) # just temporary lon values, will not affect output of this function.
    [X_map,Y_map] = np.meshgrid(lon_lim,lat_lim)
    lat_map = np.flipud(Y_map)

    return lat_map

def classifyFinalYield( est_yield):

    ''' Classifying Final Yield Map
    class 5 = very suitable = yields are equivalent to 80% or more of the overall maximum yield,
    class 4 = suitable = yields between 60% and 80%,
    class 3 = moderately suitable = yields between 40% and 60%,
    class 2 = marginally suitable = yields between 20% and 40%,
    class 1 = not suitable = yields between 0% and 20%.
    '''

    est_yield_max = np.amax( est_yield[est_yield>0] )
    est_yield_min = np.amin( est_yield[est_yield>0] )

    est_yield_20P = (est_yield_max-est_yield_min)*(20/100) + est_yield_min
    est_yield_40P = (est_yield_max-est_yield_min)*(40/100) + est_yield_min
    est_yield_60P = (est_yield_max-est_yield_min)*(60/100) + est_yield_min
    est_yield_80P = (est_yield_max-est_yield_min)*(80/100) + est_yield_min

    est_yield_class = np.zeros(est_yield.shape)

    est_yield_class[ np.all([0<est_yield, est_yield<=est_yield_20P], axis=0) ] = 1 # not suitable
    est_yield_class[ np.all([est_yield_20P<est_yield, est_yield<=est_yield_40P], axis=0) ] = 2 # marginally suitable
    est_yield_class[ np.all([est_yield_40P<est_yield, est_yield<=est_yield_60P], axis=0) ] = 3 # moderately suitable
    est_yield_class[ np.all([est_yield_60P<est_yield, est_yield<=est_yield_80P], axis=0) ] = 4 # suitable
    est_yield_class[ np.all([est_yield_80P<est_yield], axis=0)] = 5 # very suitable

    return est_yield_class

def saveRaster( ref_raster_path, out_path, numpy_raster):
    """Save NumPy arrays/matrices to GeoTIFF files

    Args:
        ref_raster_path (string): File path to referece GeoTIFF for geo-tagged info.
        out_path (string): Path for the created GeoTIFF to be saved as/to
        numpy_raster (2D NumPy): the arrays to be saveda as GeoTIFF
    """        
    # Read random image to get projection data
    img = gdal.Open(ref_raster_path)
    # allocating space in hard drive
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(out_path, img.RasterXSize, img.RasterYSize, 1, gdal.GDT_Float32)
    # set image paramenters (imfrormation related to cordinates)
    outdata.SetGeoTransform(img.GetGeoTransform())
    outdata.SetProjection(img.GetProjection())
    # write numpy matrix as new band and set no data value for the band
    outdata.GetRasterBand(1).WriteArray(numpy_raster)
    outdata.GetRasterBand(1).SetNoDataValue(-999)
    # flush data from memory to hard drive
    outdata.FlushCache()
    outdata=None

def averageRasters( raster_3d):
    """Averaging a list of raster files in time dimension

    Args:
        raster_3d (3D NumPy array): any climate data

    Returns:
        2D NumPy: the averaged climate data into 'one year' array
    """        
    # input should be a 3D raster and averaging will be done through last dimension (usually corresponding to years)
    return np.sum(raster_3d, axis=2)/raster_3d.shape[-1]

def windSpeedAt2m( wind_speed, altitude):
    """Convert windspeed at any altitude to those at 2m altitude

    Args:
        wind_speed (1D,2D,or 3D NumPy array): wind speed
        altitude (float): altitude of wind speed measurement above ground [m]

    Returns:
        1D,2D,or 3D NumPy array: Converted wind speed at 2m altitude
    """        
    # this function converts wind speed from a particular altitude of measurement to 2m altitude. wind_speed can be a numpy array (can be 1D, 2D or 3D)
    return wind_speed * (4.868/np.log(67.75*altitude-5.42))

def getYieldGap( potential, actual):
    """
    Calculates the yield gap production (Difference between
    potential yield and the actual yield production)
    
    Args:
        potential (2-D NumPy Arry, kg/ha): potential yield production
        actual (2-D NumPy Array, kg/ha): actual yield production
    """

    return np.subtract(potential, actual)
