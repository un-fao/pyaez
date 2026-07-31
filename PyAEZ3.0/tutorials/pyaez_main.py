import os
import sys
import numpy as np
from osgeo import gdal
import pandas as pd
import xarray as xr
from types import MethodType
import matplotlib.pyplot as plt
from typing import Optional, Union

# Ensure GDAL exceptions are raised
gdal.UseExceptions()

def saveRaster(
        ref_raster_path: str,
        out_path: str,
        numpy_raster: np.ndarray,
        nodata: Optional[Union[int, float]] = -999,
        compress: Optional[str] = "LZW"
    ) -> None:
    """
    Save a NumPy array to a GeoTIFF using georeferencing from a reference raster.

    This function:
      - Opens a reference GeoTIFF and copies its GeoTransform and Projection,
      - Validates that the NumPy array shape matches the reference raster size,
      - Creates (or overwrites) a single‑band GeoTIFF at `out_path`,
      - Writes the array as float or integer depending on the input dtype,
      - Sets NoData on the output band (default: -999), converting NaNs to NoData,
      - Optionally applies lossless compression (default: LZW).

    Parameters
    ----------
    ref_raster_path : str
        Path to a georeferenced reference raster (GeoTIFF). Its spatial metadata
        (GeoTransform and Projection) will be used for the output.
    out_path : str
        Destination path for the output GeoTIFF. The parent folder is created
        if it does not exist.
    numpy_raster : numpy.ndarray
        2D array of shape (rows, cols) to be saved as a single‑band GeoTIFF.
        The shape must match the reference raster’s (RasterYSize, RasterXSize).
    nodata : int or float, optional
        NoData value to set on the output band. If the array contains NaNs,
        they are replaced with this value before writing. Default: -999.
    compress : str or None, optional
        GDAL creation option for compression (e.g., "LZW", "DEFLATE").
        Use None to disable compression. Default: "LZW".

    Returns
    -------
    None
        Writes the GeoTIFF to disk; raises exceptions on failure.

    Raises
    ------
    FileNotFoundError
        If `ref_raster_path` does not exist or the raster cannot be opened.
    ValueError
        If `numpy_raster` is not 2D or its shape does not match the reference.
    OSError
        If GDAL fails to create or write the output dataset.

    Notes
    -----
    - The output data type is inferred from `numpy_raster.dtype`:
        * floating dtype → `GDT_Float32`
        * integer dtype  → `GDT_Int32` (or `GDT_Int16` for smaller ranges; customize if needed)
    - If you need to write **multiple bands**, extend this function by creating
      `driver.Create(..., bands=n_bands, ...)` and looping over bands.

    Examples
    --------
    >>> arr = np.random.random((100, 200)).astype(np.float32)
    >>> saveRaster("ref.tif", "out.tif", arr, nodata=-999, compress="LZW")
    """
    # --- Validate reference raster ---
    if not os.path.exists(ref_raster_path):
        raise FileNotFoundError(f"Reference raster not found: {ref_raster_path}")

    ref_ds = gdal.Open(ref_raster_path)
    if ref_ds is None:
        raise FileNotFoundError(f"GDAL could not open reference raster: {ref_raster_path}")

    ref_cols = ref_ds.RasterXSize
    ref_rows = ref_ds.RasterYSize
    geo_transform = ref_ds.GetGeoTransform()
    projection = ref_ds.GetProjection()

    # --- Validate array ---
    if numpy_raster.ndim != 2:
        raise ValueError(f"`numpy_raster` must be 2D (rows, cols); got shape {numpy_raster.shape}")

    rows, cols = numpy_raster.shape
    if (rows != ref_rows) or (cols != ref_cols):
        raise ValueError(
            "Array shape does not match reference raster size: "
            f"array={numpy_raster.shape} vs reference={(ref_rows, ref_cols)}"
        )

    # --- Prepare output folder ---
    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # --- Decide GDAL data type from NumPy dtype ---
    if np.issubdtype(numpy_raster.dtype, np.floating):
        gdal_dtype = gdal.GDT_Float32
        arr_to_write = numpy_raster.astype(np.float32, copy=False)
    elif np.issubdtype(numpy_raster.dtype, np.signedinteger):
        # You can choose Int16 for smaller ranges; using Int32 by default
        gdal_dtype = gdal.GDT_Int32
        arr_to_write = numpy_raster.astype(np.int32, copy=False)
    elif np.issubdtype(numpy_raster.dtype, np.unsignedinteger):
        gdal_dtype = gdal.GDT_UInt32
        arr_to_write = numpy_raster.astype(np.uint32, copy=False)
    else:
        # Fallback to float32
        gdal_dtype = gdal.GDT_Float32
        arr_to_write = numpy_raster.astype(np.float32, copy=False)

    # --- Replace NaNs with NoData if needed ---
    if np.issubdtype(arr_to_write.dtype, np.floating):
        if np.isnan(arr_to_write).any():
            arr_to_write = np.where(np.isnan(arr_to_write), nodata, arr_to_write)

    # --- Create dataset with compression option if requested ---
    driver = gdal.GetDriverByName("GTiff")
    options = []
    if compress:
        options.append(f"COMPRESS={compress}")
        # You can also add: options += ["TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"]

    out_ds = driver.Create(out_path, ref_cols, ref_rows, 1, gdal_dtype, options=options)
    if out_ds is None:
        raise OSError(f"GDAL failed to create output dataset: {out_path}")

    # --- Set spatial metadata ---
    if geo_transform:
        out_ds.SetGeoTransform(geo_transform)
    if projection:
        out_ds.SetProjection(projection)

    # --- Write data ---
    band = out_ds.GetRasterBand(1)
    if band is None:
        out_ds = None
        raise OSError("Failed to access output band 1.")
    band.WriteArray(arr_to_write)

    # --- Set NoData & flush ---
    if nodata is not None:
        band.SetNoDataValue(float(nodata))
    band.FlushCache()
    out_ds.FlushCache()

    # --- Clean up ---
    band = None
    out_ds = None

    

def initialize_clim(
        year: str,
        country_name: str,
        country_mask_name: str,
        elevation_filename: str,
        daily: bool = True,
        mask_value: int = 1,
    ):
    """
    Initialize and configure a ClimateRegime object for AEZ simulation.

    This function:
      - Validates and switches the working directory to the project root,
      - Loads climate input arrays (tasmin, tasmax, pr, hurs, sfcWind, rsds, time_array),
      - Loads the study-area mask and elevation raster,
      - Extracts geospatial bounds from the mask,
      - Filters climate data for the specified year if `daily=True`,
      - Configures the ClimateRegime object with climate and soil water data.

    Parameters
    ----------
    year : str
        Year to filter time series (e.g., "1990"). Used only if `daily=True`.
    country_name : str
        Country name (currently not used for file paths but kept for consistency).
    country_mask_name : str
        Filename of the country/study-area mask raster located in `data_input/`.
    elevation_filename : str
        Filename of the elevation raster located in `data_input/`.
    daily : bool, optional
        If True, filter climate arrays for the given year and pass daily data.
        If False, pass full arrays (monthly or climatology).
    mask_value : int, optional
        Value in the mask raster representing the study area (default: 1).

    Returns
    -------
    clim_reg : pyaez.ClimateRegime.ClimateRegime
        A configured ClimateRegime instance with:
        - study area mask set,
        - location and terrain data set,
        - daily or monthly climate and soil water data loaded.

    Notes
    -----
    - **Side effects**: Changes current working directory to the project root and appends it to `sys.path`.
    - **File dependencies**:
        * Climate arrays in `data_input/CAVA_data`:
          `tasmax.npy`, `tasmin.npy`, `pr.npy`, `hurs.npy`, `sfcWind.npy`, `rsds.npy`, `time_array.npy`.
        * Mask and elevation rasters in `data_input/`:
          `country_mask_name`, `elevation_filename`.
    - **Array shapes**: Climate arrays are expected as 3-D `(rows, cols, time)`.
      Mask and elevation rasters must align spatially with climate grids.

    Raises
    ------
    FileNotFoundError
        If the working directory, mask, elevation, or any climate input file is missing.
    ValueError
        If `daily=True` but no time indices match the specified `year`.

    Examples
    --------
       >>> clim_reg = initialize_clim(
    ...     year="1990",
    ...     country_name="Ghana",
    ...     country_mask_name="ghana_mask.tif",
    ...     elevation_filename="srtm_ghana.tif",
    ...     daily=True,
    ...     mask_value=1
    ... )
    >>> # ClimateRegime object is now ready for NB1 computations.
    """

    cwd=os.getcwd()
    work_dir = os.path.dirname(cwd)

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



def initialize_aez(
        year: str,
        country_name: str,
        country_mask_name: str,
        elevation_filename: str,
        management: str = 'HI',
        daily: bool = True
    ):
    """
    Initialize and configure an AEZ ``CropSimulation`` object for a given country/year.

    This routine:
      1) switches the current working directory to ``work_dir`` and adds it to ``sys.path``,
      2) loads climate inputs (tasmin/tasmax/pr/hurs/sfcWind/rsds and time array),
      3) loads the study-area mask and elevation raster,
      4) filters the climate series to the requested ``year`` (removing Feb 29 for leap years),
      5) configures the AEZ model with daily or monthly climate & water data,
      6) loads crop/crop-cycle parameters based on ``management``,
      7) imports NB1 screening rasters (Thermal Climate & Permafrost),
      8) sets up crop-specific temperature profile rules, and
      9) imports LGP and LGPT rasters (NB1 outputs) into the model.

    Parameters
    ----------
    year : str
        Year selector used to filter the time dimension of the climate arrays.
        Example: ``"1990"``. If the selected year has 366 daily entries, day index 59
        (Feb 29) is removed to keep 365 days.
    country_name : str
        Country (or study area) identifier used to build NB1 filenames. The function
        expects NB1 rasters with names like:
        ``{country_name.lower()}_ThermalClimate_{year}.tif`` and
        ``{country_name.lower()}_permafrost_{year}.tif`` under ``data_output/NB1/``.
    country_mask_name : str
        Filename of the country/study-area mask raster located in
        ``{work_dir}/data_input/``. Example: ``"ghana_mask.tif"``.
    elevation_filename : str
        Filename of the elevation raster located in ``{work_dir}/data_input/``.
        Must be coregistered with the mask (same extent/resolution/CRS).
    management : {'HI','LI'}, optional
        Crop parameter set. ``'HI'`` (high input) is supported and loads
        ``./data_input/input_crop_TSUM_parameters_HI.xlsx``.
        ``'LI'`` currently raises ``RuntimeError`` because low‑input parameters
        are not provided.
    daily : bool, optional
        If ``True``, pass **daily** climate series for the selected year to the model.
        If ``False``, pass **monthly** climatologies (full arrays without per‑year filtering).

    Returns
    -------
    aez : pyaez.CropSimulation.CropSimulation
        A configured AEZ CropSimulation instance with:
        - study area mask and location/terrain set,
        - daily or monthly climate & water data loaded,
        - crop & crop-cycle parameters loaded,
        - thermal climate and permafrost screening applied,
        - crop-specific temperature profile rules set,
        - LGP and LGPT rasters imported.
    
    Notes
    -----
    - **Side effects**: Changes process CWD (``os.chdir(work_dir)``) and modifies ``sys.path``.
    - **File dependencies**:
        * Climate inputs in ``{work_dir}/data_input/CAVA_data``:
          ``tasmax.npy``, ``tasmin.npy``, ``pr.npy``, ``hurs.npy``,
          ``sfcWind.npy``, ``rsds.npy``, ``time_array.npy``.
        * Mask and elevation in ``{work_dir}/data_input``:
          ``country_mask_name``, ``elevation_filename``.
        * Crop parameter workbook:
          ``./data_input/input_crop_TSUM_parameters_HI.xlsx`` (when ``management='HI'``).
        * NB1 rasters in ``{work_dir}/data_output/NB1``:
          ``{country_name.lower()}_ThermalClimate_{year}.tif``,
          ``{country_name.lower()}_permafrost_{year}.tif``,
          ``{country_name}_LGP_{year}.tif``,
          ``{country_name}_LGPt5_{year}.tif``,
          ``{country_name}_LGPt10_{year}.tif``.
    - **Array shapes**: This function assumes 3‑D arrays ``(rows, cols, time)`` for climate inputs.
      The country mask and elevation rasters must align spatially with the climate grids.

    Raises
    ------
    RuntimeError
        If ``management='LI'`` (low‑input parameters not yet available).
    NameError
        If ``management`` is neither ``'HI'`` nor ``'LI'``.
    FileNotFoundError / OSError
        If any required input file cannot be opened/read by GDAL or NumPy.
    ValueError
        If time filtering yields an unexpected number of daily slices.

    Examples
    --------
    >>> aez = initialize_aez(
    ...     year="1990",
    ...     country_name="Ghana",
    ...     country_mask_name="ghana_mask.tif",
    ...     elevation_filename="srtm_ghana.tif",
    ...     management="HI",
    ...     daily=True,
    ... )
    >>> # The model is now configured and ready for simulation steps.
    """

    # Set working directory
    cwd=os.getcwd()
    work_dir = os.path.dirname(cwd)
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

    if management == 'HI':
        crop_param_filename = r'../data_input/input_crop_TSUM_parameters_HI.xlsx'
    elif management == 'LI':
        # crop_param_filename = r'./data_input/input_crop_TSUM_parameters_LI.xlsx'
        raise RuntimeError("Low input parameters are yet to be provided")
    else:
        raise NameError("The variable management can be either 'HI' or 'LI'.")

    df = pd.read_excel(crop_param_filename, sheet_name='parameters')
    
    # Filter non-NaN crop names
    available_crops = df[df['Crop_name'].notna()]['Crop_name'].unique()
    
    print("Available crops:")
    for i, crop in enumerate(available_crops, start=1):
        print(f"{i}. {crop}")
    
    choice = int(input("Select a crop by number: "))
    crop_name = available_crops[choice - 1]
    print(f"You selected: {crop_name}")

    aez.readCropandCropCycleParameters(file_path = crop_param_filename, 
                                   crop_name = crop_name)

    tclimate = gdal.Open(rf"../data_output/NB1/{country_name.lower()}_ThermalClimate_{year}.tif").ReadAsArray()
    permafrost_class = gdal.Open(rf"../data_output/NB1/{country_name.lower()}_permafrost_{year}.tif").ReadAsArray()
    
    
    # Thermal Climate screening
    # aez.setThermalClimateScreening(tclimate, no_t_climate=[6,7,8,9,10,11,12])
    aez.setThermalClimateScreening(tclimate, no_t_climate=[])
    
    # New Thermal Screening: Permafrost Screening
    aez.setPermafrostScreening(permafrost_class= permafrost_class)
    
    # Updated Temperature Profile screenign routine
    aez.setupCropSpecificRule(file_path = r'../data_input/crop-specific_rule.xlsx',
                               crop_name = crop_name.split('_')[0])

    lgp = gdal.Open(rf'../data_output/NB1/{country_name}_LGP_{year}.tif').ReadAsArray()
    lgpt0 = gdal.Open(rf'../data_output/NB1/{country_name}_LGPt0_{year}.tif').ReadAsArray()
    lgpt5 = gdal.Open(rf'../data_output/NB1/{country_name}_LGPt5_{year}.tif').ReadAsArray()
    lgpt10 = gdal.Open(rf'../data_output/NB1/{country_name}_LGPt10_{year}.tif').ReadAsArray()
    lgp_equv = gdal.Open(rf'../data_output/NB1/{country_name}_LGPEquivalent_{year}.tif').ReadAsArray()
    
    
    aez.ImportLGPandLGPT(lgp = lgp, lgpt0 = lgpt0, lgpt5 = lgpt5, lgpt10= lgpt10)

    return aez



def classifyFinalYield(est_yield):
    """
    Classify estimated crop yield values into suitability classes.

    This function assigns each cell in the estimated yield map to one of five
    suitability classes based on its relative position between the minimum and
    maximum positive yield values:

        Class 5 (Very suitable): yield >= 80% of the range above min yield
        Class 4 (Suitable):      yield >= 60% and < 80%
        Class 3 (Moderate):      yield >= 40% and < 60%
        Class 2 (Marginal):      yield >= 20% and < 40%
        Class 1 (Not suitable):  yield > 0% and < 20%

    Parameters
    ----------
    est_yield : numpy.ndarray
        2D array of estimated yields (numeric). Zero or negative values
        are ignored when computing thresholds.

    Returns
    -------
    est_yield_class : numpy.ndarray
        Array of the same shape as `est_yield` with integer class codes:
        {0: no data or non-positive yield, 1–5: suitability classes}.

    Notes
    -----
    - Thresholds are computed using the min and max of positive yields only.
    - Classification uses inclusive upper bounds for each range.
    - Output classes:
        0 = no data / yield <= 0
        1 = not suitable
        2 = marginally suitable
        3 = moderately suitable
        4 = suitable
        5 = very suitable

    Example
    -------
    >>> import numpy as np
    >>> yields = np.array([[0, 1.5, 3.0],
    ...                    [4.5, 6.0, 8.0]])
    >>> classifyFinalYield(yields)
    array([[0., 1., 2.],
           [3., 4., 5.]])
    """

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


def initialize_climatic_constraints(
        work_dir: str,
        country_name: str,
        country_mask_name: str,
        elevation_filename: str,
        management: str,
        year: str
    ):
    """
    Initialize and configure a ClimaticConstraints object for AEZ simulation.

    This function:
      - Changes the working directory to `work_dir` and adds it to `sys.path`,
      - Loads the study-area mask and elevation raster,
      - Extracts geospatial bounds from the mask,
      - Loads crop parameter workbook based on `management` (HI or LI),
      - Prompts the user to select a crop interactively,
      - Loads climate input arrays (tasmin, tasmax, pr, hurs, sfcWind, rsds, time_array),
      - Filters climate data for the specified `year` (removing Feb 29 for leap years),
      - Initializes and configures the ClimaticConstraints object with mask, terrain, and climate data,
      - Attaches the selected crop name to the object.

    Parameters
    ----------
    work_dir : str
        Absolute or relative path to the project working directory.
        Expected structure includes:
        - `data_input/` (mask, elevation, crop parameters, climate arrays),
        - `data_output/NB3/` (output folder created if missing).
    country_name : str
        Country or study-area name (used for file naming and logging).
    country_mask_name : str
        Filename of the country/study-area mask raster located in `data_input/`.
    elevation_filename : str
        Filename of the elevation raster located in `data_input/`.
    management : {'HI','LI'}
        Crop parameter set:
        - 'HI' loads `input_crop_TSUM_parameters_HI.xlsx`,
        - 'LI' currently raises `RuntimeError` (not implemented).
    year : str
        Year selector used to filter the time dimension of climate arrays.
        Example: `"1990"`. If the selected year has 366 daily entries, day index 59
        (Feb 29) is removed to keep 365 days.

    Returns
    -------
    clim_con : pyaez.ClimaticConstraints.ClimaticConstraints
        A configured ClimaticConstraints instance with:
        - study area mask and terrain set,
        - daily climate data for the selected year loaded,
        - crop name attached.

    Notes
    -----
    - **Side effects**:
        * Changes current working directory to `work_dir`,
        * Appends `work_dir` to `sys.path`,
        * Creates `data_output/NB3` if missing,
        * Prompts user for crop selection (interactive).
    - **File dependencies**:
        * Mask and elevation in `data_input/`:
          `country_mask_name`, `elevation_filename`.
        * Crop parameter workbook:
          `input_crop_TSUM_parameters_HI.xlsx` (when `management='HI'`).
        * Climate arrays in `data_input/CAVA_data`:
          `tasmax.npy`, `tasmin.npy`, `pr.npy`, `hurs.npy`,
          `sfcWind.npy`, `rsds.npy`, `time_array.npy`.
    - **Array shapes**:
        Climate arrays are expected as 3-D `(rows, cols, time)`.
        Mask and elevation rasters must align spatially with climate grids.

    Raises
    ------
    FileNotFoundError
        If the working directory, mask, elevation, or any climate input file is missing.
    RuntimeError
        If `management='LI'` (low-input parameters not yet available).
    NameError
        If `management` is neither 'HI' nor 'LI'.
    ValueError
        If time filtering yields an unexpected number of daily slices.
    IndexError
        If the crop selection index is invalid.

    Examples
    --------
    >>> clim_con = initialize_climatic_constraints(
    ...     work_dir=r"D:/Projects/AEZ",
    ...     country_name="Ghana",
    ...     country_mask_name="ghana_mask.tif",
    ...     elevation_filename="srtm_ghana.tif",
    ...        ...     management="HI",
    ...     year="1990"
    ... )
    >>> # ClimaticConstraints object is now ready for NB3 computations.
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
    if management == 'HI':
        crop_param_filename = r'./data_input/input_crop_TSUM_parameters_HI.xlsx'
    elif management == 'LI':
        # crop_param_filename = r'./data_input/input_crop_TSUM_parameters_LI.xlsx'
        raise RuntimeError("Low input parameters are yet to be provided")
    else:
        raise NameError("The variable management can be either 'HI' or 'LI'.")
        
    df = pd.read_excel(crop_param_filename, sheet_name='parameters')
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



def compute_yield(
        clim_con,
        condition_type: str,
        work_dir: str,
        country_name: str,
        year,
        management: str = 'HI',
        plot_results: bool = False
    ):
    """
    Compute climate-adjusted yield and reduction factor for rainfed or irrigated
    conditions, selecting the correct parameter file based on the crop name and
    management level (HI/LI).

    This function:
      - Loads the baseline yield raster (NB2) for the requested condition,
      - Loads agro-climatic indicators (LGP, LGPT10, LGPEquivalent) from NB1,
      - Selects the appropriate parameter workbook (Module 3) by `crop_name` and `management`,
      - Applies climatic reduction factors and constraints via `clim_con`,
      - Returns the climate-adjusted yield and the overall reduction factor (Fc3),
      - Optionally plots original yield, constrained yield, and Fc3.

    Parameters
    ----------
    clim_con : pyaez.ClimaticConstraints.ClimaticConstraints
        Configured climatic constraints object. Must provide:
        `crop_name` attribute and methods:
        `setReductionFactors(file_path)`,
        `applyClimaticConstraints(...)`,
        `getClimateAdjustedYield()`,
        `getClimateReductionFactor()`.
    condition_type : {'rainfed', 'irrigated'}
        Production condition to compute. Determines which baseline yield raster
        is loaded (NB2: `yld_rain` vs `yld_irr`) and which parameter file is used.
    work_dir : str
        Project working directory containing `data_input/` and `data_output/`.
    country_name : str
        Country/study-area identifier used in NB1/NB2 filenames.
    year : int or str
        Year used in NB1/NB2 filenames. Accepts either `int` (e.g., 1990) or `str` (e.g., "1990").
    management : {'HI','LI'}, optional
        Parameter set selection:
        - 'HI' → High-input workbooks (Module 3),
        - 'LI' → Low-input workbooks (Module 3).
        Raises `NameError` if not one of {'HI','LI'}.
    plot_results : bool, optional
        If True, displays three plots (original yield, climate-adjusted yield, Fc3).

    Returns
    -------
    clim_yield : numpy.ndarray
        Climate-constrained yield map (same shape as input yield).
    fc3 : numpy.ndarray
        Overall climatic reduction factor (Fc3), values typically in [0, 1].

    Notes
    -----
    - **File dependencies**:
        * NB2 baseline yield rasters in `data_output/NB2`:
          `{country_name}_{crop_name}_yld_rain_{year}.tif` or `{country_name}_{crop_name}_yld_irr_{year}.tif`
        * NB1 agro-climatic indicators in `data_output/NB1`:
          `{country_name}_LGP_{year}.tif`,
          `{country_name}_LGPt10_{year}.tif`,
          `{country_name}_LGPEquivalent_{year}.tif`
        * Module 3 parameter workbooks in `data_input/input_module3/`:
          crop-specific `.xlsx` files (HI or LI set)
    - **Crop name normalization**: The function derives `crop_key` from `clim_con.crop_name`
      by dropping trailing `_H` or `_L` and lowercasing (e.g., `"Maize_HI"` → `"maize"`).
    - **Array shapes**: All rasters are expected to be co-registered (same extent/resolution/CRS).
    - **Plotting**: Uses Matplotlib; `vmax` is set to the maximum of original vs. constrained yield.

    Raises
    ------
    ValueError
        If `condition_type` is not one of {'rainfed','irrigated'} or if `crop_key` is unsupported.
    NameError
        If `management` is neither 'HI' nor 'LI'.
    FileNotFoundError / OSError
        If any of the required rasters or parameter files cannot be opened.

    Examples
    --------
    >>> clim_yield, fc3 = compute_yield(
    ...     clim_con=clim_con,
    ...     condition_type='rainfed',
    ...     work_dir=r'D:/Projects/AEZ',
    ...     country_name='Ghana',
    ...     year=1990,
    ...     management='HI',
    ...     plot_results=True
    ... )
    >>> clim_yield.shape, fc3.min(), fc3.max()
       ((rows, cols), 0.0, 1.0)
    """
    
    crop_name = clim_con.crop_name

    # Validate condition type
    if condition_type not in ["rainfed", "irrigated"]:
        raise ValueError("condition_type must be 'rainfed' or 'irrigated'")

    # Map crop names to file paths
    # Map crop names to file paths
    if management == 'HI':
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
    elif management == 'LI':
        crop_files = {
            'maize': (
                r'./data_input/input_module3/Maize (tropical lowland  cultivars)_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Maize (tropical lowland  cultivars)_Low_rf_lst.xlsx'
            ),
            'rice': (
                r'./data_input/input_module3/Indica wetland rice_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Indica wetland rice_Low_rf_lst.xlsx'
            ),
            'soybean': (
                r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Soybean (tropical and subtropical cultivars)_Low_rf_lst.xlsx'
            ),
            'cassava': (
                r'./data_input/input_module3/Cassava_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Cassava_Low_rf_lst.xlsx'
            ),
            'cashew': (
                r'./data_input/input_module3/Cashew_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Cashew_Low_rf_lst.xlsx'
            ),
            'cocoa': (
                r'./data_input/input_module3/Cocoa_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Cocoa_Low_rf_lst.xlsx'
            ),
            'coffee': (
                r'./data_input/input_module3/Coffee_robusta_Low_ir_lst.xlsx',
                r'./data_input/input_module3/Coffee_robusta_Low_rf_lst.xlsx'
            )
        }
    else:
        raise NameError("The variable management can be either 'HI' or 'LI'.")

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
    crop_key = clim_con.crop_name.replace('_L', '').lower()

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
