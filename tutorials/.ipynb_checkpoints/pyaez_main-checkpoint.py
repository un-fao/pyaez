import os
import sys
import numpy as np
from osgeo import gdal
import pandas as pd
import xarray as xr
from types import MethodType

# Ensure GDAL exceptions are raised
gdal.UseExceptions()
os.environ['PROJ_LIB'] = '/opt/conda/share/proj'

def initialize_clim(
    work_dir: str,
    year: str,
    country_name: str,
    country_mask_name: str,
    elevation_filename: str,
    daily: bool = True,
    mask_value: int = 1,
):
    """
    Initialize and configure the ClimateRegime object for AEZ simulation.

    Parameters:
        work_dir (str): Working directory path.
        year (str, optional): Year to filter time series (currently unused).
        country_name (str, optional): Country name (currently unused).
        country_mask_name (str): Filename of the country mask raster.
        elevation_filename (str): Filename of the elevation raster.
        daily (bool): Whether to use daily or monthly data.
        mask_value (int): Value in mask representing the study area.

    Returns:
        ClimateRegime: Configured ClimateRegime object.
    """
    # Validate working directory
    if not os.path.isdir(work_dir):
        raise FileNotFoundError(f"Working directory not found: {work_dir}")

    os.chdir(work_dir)
    sys.path.append(work_dir)

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB1')
    os.makedirs(folder_path, exist_ok=True)

    from pyaez import ClimateRegime

    # Initialize ClimateRegime object
    clim_reg = ClimateRegime.ClimateRegime()

    # Load climate data
    data_input = os.path.join(work_dir, 'data_input', 'CAVA_data')
    climate_files = {
        "tasmax": "tasmax.npy",
        "tasmin": "tasmin.npy",
        "pr": "pr.npy",
        "hurs": "hurs.npy",
        "sfcWind": "sfcWind.npy",
        "rsds": "rsds.npy",
        "time_array": "time_array.npy"
    }

    # Check all files exist
    for key, fname in climate_files.items():
        path = os.path.join(data_input, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing climate data file: {path}")

    max_temp = np.load(os.path.join(data_input, climate_files["tasmax"]))
    min_temp = np.load(os.path.join(data_input, climate_files["tasmin"]))
    precipitation = np.load(os.path.join(data_input, climate_files["pr"]))
    rel_humidity = np.load(os.path.join(data_input, climate_files["hurs"]))
    wind_speed = np.load(os.path.join(data_input, climate_files["sfcWind"]))
    short_rad = np.load(os.path.join(data_input, climate_files["rsds"]))
    times = np.load(os.path.join(data_input, climate_files["time_array"]))

    # Load mask and elevation
    mask_path = os.path.join(work_dir, 'data_input', country_mask_name)
    if not os.path.exists(mask_path):
        raise FileNotFoundError(f"Mask file not found: {mask_path}")

    mask_gdal = gdal.Open(mask_path)
    mask = mask_gdal.ReadAsArray()

    elevation_path = os.path.join(work_dir, 'data_input', elevation_filename)
    if not os.path.exists(elevation_path):
        raise FileNotFoundError(f"Elevation file not found: {elevation_path}")

    elevation = gdal.Open(elevation_path).ReadAsArray()

    # Extract geotransform and bounds
    geotransform = mask_gdal.GetGeoTransform()
    minx = geotransform[0]
    miny = geotransform[3] + geotransform[5] * mask_gdal.RasterYSize
    maxx = geotransform[0] + geotransform[1] * mask_gdal.RasterXSize
    maxy = geotransform[3]

    lat_min, lat_max = miny, maxy

    # Configure ClimateRegime
    clim_reg.setStudyAreaMask(mask, mask_value)
    clim_reg.setLocationTerrainData(lat_min, lat_max, elevation)

    red_times_idx = [i for i,t in enumerate(times) if year in str(t)]
    
    # Set climate and soil water data
    if daily:
        if red_times_idx is None:
            raise ValueError("red_times_idx must be provided for daily data.")
        clim_reg.setClimateAndSoilWaterData(
            min_temp[:, :, red_times_idx],
            max_temp[:, :, red_times_idx],
            precipitation[:, :, red_times_idx],
            short_rad[:, :, red_times_idx],
            wind_speed[:, :, red_times_idx],
            rel_humidity[:, :, red_times_idx],
            itflg=5
        )
    else:
        clim_reg.setClimateAndSoilWaterData(
            min_temp, max_temp, precipitation, short_rad, wind_speed, rel_humidity
        )

    return clim_reg


def initialize_aez(work_dir, year, country_name, country_mask_name, elevation_filename, daily=True):
    """
    Initialize and configure the AEZ CropSimulation object.
    
    Parameters:
        work_dir (str): Working directory path.
        year (str): Year to filter time series.
        daily (bool): Whether to use daily or monthly data.
    
    Returns:
        aez (CropSimulation.CropSimulation): Configured AEZ object.
    """
    # Set working directory
    os.chdir(work_dir)
    sys.path.append(work_dir)

    # Import your library
    from pyaez import CropSimulation
    from pyaez.UtilitiesCalc import saveRaster, classifyFinalYield

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB2')
    os.makedirs(folder_path, exist_ok=True)

    # Initialize CropSimulation
    aez = CropSimulation.CropSimulation()

    # Load climate data
    data_input = os.path.join(work_dir, 'data_input', 'CAVA_data')
    max_temp = np.load(os.path.join(data_input, 'tasmax.npy'))
    min_temp = np.load(os.path.join(data_input, 'tasmin.npy'))
    precipitation = np.load(os.path.join(data_input, 'pr.npy'))
    rel_humidity = np.load(os.path.join(data_input, 'hurs.npy'))
    wind_speed = np.load(os.path.join(data_input, 'sfcWind.npy'))
    short_rad = np.load(os.path.join(data_input, 'rsds.npy'))
    times = np.load(os.path.join(data_input, 'time_array.npy'))

    # Load mask and elevation
    mask_path = os.path.join(work_dir, 'data_input', country_mask_name)
    mask_gdal = gdal.Open(mask_path)
    mask = mask_gdal.ReadAsArray()
    elevation = gdal.Open(os.path.join(work_dir, 'data_input', elevation_filename)).ReadAsArray()

    # Extract geotransform and bounds
    geotransform = mask_gdal.GetGeoTransform()
    minx = geotransform[0]
    miny = geotransform[3] + geotransform[5] * mask_gdal.RasterYSize
    maxx = geotransform[0] + geotransform[1] * mask_gdal.RasterXSize
    maxy = geotransform[3]

    lat_min, lat_max = miny, maxy

    # Filter times for the given year
    red_times_idx = [i for i, t in enumerate(times) if year in str(t)]
    if len(red_times_idx) == 366:  # Leap year adjustment
        red_times_idx = np.delete(red_times_idx, 59)

    # Configure AEZ object
    aez.setStudyAreaMask(admin_mask=mask, no_data_value=0)
    aez.setLocationTerrainData(lat_min=lat_min, lat_max=lat_max, elevation=elevation)

    # Soil water holding capacity and rooting depth
    Sa, D = 100., 1.

    if daily:
        aez.setDailyClimateAndWaterData(
            min_temp[:, :, red_times_idx],
            max_temp[:, :, red_times_idx],
            precipitation[:, :, red_times_idx],
            short_rad[:, :, red_times_idx],
            wind_speed[:, :, red_times_idx],
            rel_humidity[:, :, red_times_idx],
            Sa, D
        )
    else:
        aez.setMonthlyClimateAndWaterData(
            min_temp, max_temp, precipitation, short_rad, wind_speed, rel_humidity, Sa, D
        )

    df = pd.read_excel(r'./data_input/input_crop_TSUM_parameters_TOGO_HI.xlsx')

    # Filter non-NaN crop names
    available_crops = df[df['Crop_name'].notna()]['Crop_name'].unique()
    
    print("Available crops:")
    for i, crop in enumerate(available_crops, start=1):
        print(f"{i}. {crop}")
    
    choice = int(input("Select a crop by number: "))
    crop_name = available_crops[choice - 1]
    print(f"You selected: {crop_name}")

    aez.readCropandCropCycleParameters(file_path = r'./data_input/input_crop_TSUM_parameters_TOGO_HI.xlsx', 
                                   crop_name = crop_name)

    tclimate = gdal.Open(rf"./data_output/NB1/togo_ThermalClimate_{year}.tif").ReadAsArray()
    permafrost_class = gdal.Open(rf"./data_output/NB1/togo_permafrost_{year}.tif").ReadAsArray()
    
    
    # Thermal Climate screening
    # aez.setThermalClimateScreening(tclimate, no_t_climate=[6,7,8,9,10,11,12])
    aez.setThermalClimateScreening(tclimate, no_t_climate=[])
    
    # New Thermal Screening: Permafrost Screening
    aez.setPermafrostScreening(permafrost_class= permafrost_class)
    
    # Updated Temperature Profile screenign routine
    aez.setupCropSpecificRule(file_path = r'./data_input/crop-specific_rule_TOGO.xlsx',
                               crop_name = crop_name.split('_')[0])

    lgp = gdal.Open(rf'./data_output/NB1/{country_name}_LGP_{year}.tif').ReadAsArray()
    lgpt5 = gdal.Open(rf'./data_output/NB1/{country_name}_LGPt5_{year}.tif').ReadAsArray()
    lgpt10 = gdal.Open(rf'./data_output/NB1/{country_name}_LGPt10_{year}.tif').ReadAsArray()
    lgp_equv = gdal.Open(rf'./data_output/NB1/{country_name}_LGPEquivalent_{year}.tif').ReadAsArray()
    
    
    aez.ImportLGPandLGPT(lgp = lgp, lgpt5 = lgpt5, lgpt10= lgpt10)

    return aez


def classifyFinalYield(est_yield):

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

import os
import sys
import numpy as np
import pandas as pd
from osgeo import gdal

def initialize_climatic_constraints(work_dir, country_name, country_mask_name, elevation_filename, year):
    """
    Initialize ClimaticConstraints object and prepare all required data.

    Parameters:
        work_dir (str): Working directory path.
        country_name (str): Country name for file naming.
        country_mask_name (str): Filename of the country mask raster.
        elevation_filename (str): Filename of the elevation raster.
        year (str): Year for filtering time series.

    Returns:
        clim_con (ClimaticConstraints.ClimaticConstraints): Configured climatic constraints object.
        crop_name (str): Selected crop name.
        yield_map_rain (ndarray): Yield map for rainfed conditions.
        yield_map_irr (ndarray): Yield map for irrigated conditions.
        lgp, lgp_equv, lgp10 (ndarray): Agro-climatic indicators.
    """
    # Ensure working directory
    os.chdir(work_dir)
    if work_dir not in sys.path:
        sys.path.insert(0, work_dir)

    # Import AEZ library
    from pyaez import ClimaticConstraints

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB3')
    os.makedirs(folder_path, exist_ok=True)

    # Load mask and elevation
    mask_path = os.path.join(work_dir, 'data_input', country_mask_name)
    mask_gdal = gdal.Open(mask_path)
    mask = mask_gdal.ReadAsArray()
    elevation = gdal.Open(os.path.join(work_dir, 'data_input', elevation_filename)).ReadAsArray()

    # Extract geotransform and bounds
    geotransform = mask_gdal.GetGeoTransform()
    lat_min = geotransform[3] + geotransform[5] * mask_gdal.RasterYSize
    lat_max = geotransform[3]

    # Load crop list and let user choose
    df = pd.read_excel(os.path.join(work_dir, 'data_input', 'input_crop_TSUM_parameters_TOGO_HI.xlsx'))
    available_crops = df[df['Crop_name'].notna()]['Crop_name'].unique()

    print("Available crops:")
    for i, crop in enumerate(available_crops, start=1):
        print(f"{i}. {crop}")
    choice = int(input("Select a crop by number: "))
    crop_name = available_crops[choice - 1]
    print(f"You selected: {crop_name}")

    # Load climate data
    data_input = os.path.join(work_dir, 'data_input', 'CAVA_data')
    max_temp = np.load(os.path.join(data_input, 'tasmax.npy'))
    min_temp = np.load(os.path.join(data_input, 'tasmin.npy'))
    precipitation = np.load(os.path.join(data_input, 'pr.npy'))
    rel_humidity = np.load(os.path.join(data_input, 'hurs.npy'))
    wind_speed = np.load(os.path.join(data_input, 'sfcWind.npy'))
    short_rad = np.load(os.path.join(data_input, 'rsds.npy'))
    times = np.load(os.path.join(data_input, 'time_array.npy'))

    
    # Filter times for the given year
    red_times_idx = [i for i, t in enumerate(times) if year in str(t)]
    if len(red_times_idx) == 366:
        red_times_idx = np.delete(red_times_idx, 59)

    # Initialize ClimaticConstraints object
    clim_con = ClimaticConstraints.ClimaticConstraints(
        lat_min=lat_min, lat_max=lat_max, elevation=elevation, mask=mask, no_mask_value=0
    )

    # Set climate data
    clim_con.setClimateData(
        min_temp=min_temp[:, :, red_times_idx],
        max_temp=max_temp[:, :, red_times_idx],
        wind_speed=wind_speed[:, :, red_times_idx],
        short_rad=short_rad[:, :, red_times_idx],
        rel_humidity=rel_humidity[:, :, red_times_idx],
        precip=precipitation[:, :, red_times_idx]
    )

    
    if not hasattr(clim_con, 'crop_name'):
        setattr(clim_con, 'crop_name', crop_name)
    else:
        clim_con.crop_name = crop_name

    
    return clim_con
   

def compute_yield(clim_con, condition_type, work_dir, country_name, year):
    """
    Compute climate-adjusted yield and reduction factor for rainfed or irrigated conditions,
    automatically selecting the correct input file based on crop name.

    Parameters:
        clim_con: Climatic constraints object (with methods setReductionFactors, applyClimaticConstraints, etc.)
        condition_type (str): 'rainfed' or 'irrigated'.

    Returns:
        tuple: (clim_yield, fc3) for the specified condition.
    """
    crop_name = clim_con.crop_name

    # Validate condition type
    if condition_type not in ["rainfed", "irrigated"]:
        raise ValueError("condition_type must be 'rainfed' or 'irrigated'")

    # Map crop names to file paths
    crop_files = {
        'maize': (
            r'./data_input/input_module3/Maize (tropical lowland  cultivars)_High_ir_lst.xlsx',
            r'./data_input/input_module3/Maize (tropical lowland  cultivars)_High_rf_lst.xlsx'
        ),
        'rice': (
            r'./data_input/input_module3/Indica wetland rice_High_ir_lst.xlsx',
            r'./data_input/input_module3/Indica wetland rice_High_rf_lst.xlsx'
        ),
        'soybean': (
            r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_High_ir_lst.xlsx',
            r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_High_rf_lst.xlsx'
        ),
        'cassava': (
            r'./data_input/input_module3/Cassava_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cassava_High_rf_lst.xlsx'
        ),
        'cashew': (
            r'./data_input/input_module3/Cashew_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cashew_High_rf_lst.xlsx'
        ),
        'cocoa': (
            r'./data_input/input_module3/Cocoa_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cocoa_High_rf_lst.xlsx'
        ),
        'coffee': (
            r'./data_input/input_module3/Coffee_robusta_High_ir_lst.xlsx',
            r'./data_input/input_module3/Coffee_robusta_High_rf_lst.xlsx'
        )
    }

    # Yield maps
    if condition_type == "rainfed":
        yield_map = gdal.Open(os.path.join(work_dir, f'data_output/NB2/{country_name}_{crop_name}_yld_rain_{year}.tif')).ReadAsArray()
    elif condition_type == "irrigated":
        yield_map = gdal.Open(os.path.join(work_dir, f'data_output/NB2/{country_name}_{crop_name}_yld_irr_{year}.tif')).ReadAsArray()

    # Agro-climatic indicators
    lgp = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGP_{year}.tif')).ReadAsArray()
    lgp10 = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPt10_{year}.tif')).ReadAsArray()
    lgp_equv = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPEquivalent_{year}.tif')).ReadAsArray()


    # Normalize crop name (remove suffix like _H)
    crop_key = clim_con.crop_name.replace('_H', '').lower()

    if crop_key not in crop_files:
        raise ValueError(f"Unsupported crop name: {crop_key}")

    filename_ir, filename_rf = crop_files[crop_key]
    filename = filename_ir if condition_type == "irrigated" else filename_rf

    # Set reduction factors
    clim_con.setReductionFactors(file_path=filename)

    # Apply climatic constraints
    clim_con.applyClimaticConstraints(
        yield_input=yield_map,
        lgp=lgp,
        lgp_equv=lgp_equv,
        lgpt10=lgp10,
        omit_yld_0=True
    )

    # Get climate-adjusted yield and reduction factor
    clim_yield = clim_con.getClimateAdjustedYield()
    fc3 = clim_con.getClimateReductionFactor()

    print(f"Computed {condition_type} yield for {crop_key} using {filename}")

    '''visualize results'''
    plt.figure(1, figsize=(22,9))
    plt.subplot(1,3,1)
    plt.imshow(yield_map_irr, vmax = np.max([clim_yield_irr, yield_map_irr]))
    plt.colorbar(shrink=0.6)
    plt.title(f'Original Irrigated Yield {crop_name}')
    
    plt.subplot(1,3,2)
    plt.imshow(clim_yield_irr_maiz, vmax = np.max([clim_yield_irr, yield_map_irr]))
    plt.colorbar(shrink=0.6)
    plt.title(f'Climate Constrainted Irrigated Yield {crop_name}')
    
    plt.subplot(1,3,3)
    plt.imshow(fc3_maiz_irr, vmax = 1)
    plt.colorbar(shrink=0.6)
    plt.title(f'Fc3 {crop_name} Irrigated')

    return clim_yield, fc3

import os
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal

def compute_yield(clim_con, condition_type, work_dir, country_name, year, plot_results=False):
    """
    Compute climate-adjusted yield and reduction factor for rainfed or irrigated conditions,
    automatically selecting the correct input file based on crop name.

    Parameters:
        clim_con: Climatic constraints object (with methods setReductionFactors, applyClimaticConstraints, etc.)
        condition_type (str): 'rainfed' or 'irrigated'.
        work_dir (str): Working directory containing input/output data.
        country_name (str): Country name for file paths.
        year (int): Year for file paths.
        plot_results (bool): If True, plots the original yield, adjusted yield, and reduction factor.

    Returns:
        tuple: (clim_yield, fc3) for the specified condition.
    """
    crop_name = clim_con.crop_name

    # Validate condition type
    if condition_type not in ["rainfed", "irrigated"]:
        raise ValueError("condition_type must be 'rainfed' or 'irrigated'")

    # Map crop names to file paths
    crop_files = {
        'maize': (
            r'./data_input/input_module3/Maize (tropical lowland  cultivars)_High_ir_lst.xlsx',
            r'./data_input/input_module3/Maize (tropical lowland  cultivars)_High_rf_lst.xlsx'
        ),
        'rice': (
            r'./data_input/input_module3/Indica wetland rice_High_ir_lst.xlsx',
            r'./data_input/input_module3/Indica wetland rice_High_rf_lst.xlsx'
        ),
        'soybean': (
            r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_High_ir_lst.xlsx',
            r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_High_rf_lst.xlsx'
        ),
        'cassava': (
            r'./data_input/input_module3/Cassava_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cassava_High_rf_lst.xlsx'
        ),
        'cashew': (
            r'./data_input/input_module3/Cashew_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cashew_High_rf_lst.xlsx'
        ),
        'cocoa': (
            r'./data_input/input_module3/Cocoa_High_ir_lst.xlsx',
            r'./data_input/input_module3/Cocoa_High_rf_lst.xlsx'
        ),
        'coffee': (
            r'./data_input/input_module3/Coffee_robusta_High_ir_lst.xlsx',
            r'./data_input/input_module3/Coffee_robusta_High_rf_lst.xlsx'
        )
    }

    # Load yield map
    if condition_type == "rainfed":
        yield_map = gdal.Open(os.path.join(work_dir, f'data_output/NB2/{country_name}_{crop_name}_yld_rain_{year}.tif')).ReadAsArray()
    else:
        yield_map = gdal.Open(os.path.join(work_dir, f'data_output/NB2/{country_name}_{crop_name}_yld_irr_{year}.tif')).ReadAsArray()

    # Load agro-climatic indicators
    lgp = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGP_{year}.tif')).ReadAsArray()
    lgp10 = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPt10_{year}.tif')).ReadAsArray()
    lgp_equv = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPEquivalent_{year}.tif')).ReadAsArray()

    # Normalize crop name
    crop_key = crop_name.replace('_H', '').lower()
    if crop_key not in crop_files:
        raise ValueError(f"Unsupported crop name: {crop_key}")

    filename_ir, filename_rf = crop_files[crop_key]
    filename = filename_ir if condition_type == "irrigated" else filename_rf

    # Apply reduction factors and constraints
    clim_con.setReductionFactors(file_path=filename)
    clim_con.applyClimaticConstraints(
        yield_input=yield_map,
        lgp=lgp,
        lgp_equv=lgp_equv,
        lgpt10=lgp10,
        omit_yld_0=True
    )

    # Get results
    clim_yield = clim_con.getClimateAdjustedYield()
    fc3 = clim_con.getClimateReductionFactor()

    print(f"Computed {condition_type} yield for {crop_key} using {filename}")

    # Optional plotting
    if plot_results:
        vmax_yield = np.max([np.nanmax(clim_yield), np.nanmax(yield_map)])
        plt.figure(figsize=(22, 9))

        plt.subplot(1, 3, 1)
        plt.imshow(yield_map, vmax=vmax_yield)
        plt.colorbar(shrink=0.6)
        plt.title(f'Original {condition_type.capitalize()} Yield ({crop_name})')

        plt.subplot(1, 3, 2)
        plt.imshow(clim_yield, vmax=vmax_yield)
        plt.colorbar(shrink=0.6)
        plt.title(f'Climate-Constrained Yield ({crop_name})')

        plt.subplot(1, 3, 3)
        plt.imshow(fc3, vmax=1)
        plt.colorbar(shrink=0.6)
        plt.title(f'Reduction Factor Fc3 ({crop_name})')

        plt.tight_layout()
        plt.show()

    return clim_yield, fc3
