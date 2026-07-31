"""
PyAEZ version 4.0 (6 June 2025)
This function library provides support for the execution of Module 0
related to the preparation of input data for Module 1.
6 June 2025: Dario Spiller, initialization of beta version

Modifications:
"""

import pandas as pd
import geopandas as gpd
from Levenshtein import distance as lev_distance  
import os
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling
import matplotlib.pyplot as plt

def get_UN_country(country_name, thr=3, return_gdf=False):
    """
    Retrieve the geometry of a country from a UN shapefile using exact or fuzzy name matching.

    Parameters:
    -----------
    country_name : str
        The name of the country to search for.
    thr : int, optional (default=3)
        The maximum Levenshtein distance allowed for fuzzy matching if no exact match is found.
    return_gdf : bool, optional (default=False)
        If True, returns both the matched country geometry and the full GeoDataFrame.

    Returns:
    --------
    AOI : GeoDataFrame or None
        A GeoDataFrame containing the dissolved geometry of the matched country, or None if no match is found.
    gdf : GeoDataFrame (optional)
        The full GeoDataFrame from the shapefile, returned only if `return_gdf` is True.
    """
    
    # Path to the UN shapefile
    UN_maps = r'data_input/UN_world_shapefile/BNDA_A1.shp'
    
    # Read the shapefile into a GeoDataFrame
    gdf = gpd.read_file(UN_maps)

    # Define a helper function to compute Levenshtein distance
    def get_distance(row_name):
        try:
            return lev_distance(row_name, country_name)
        except Exception:
            return None

    # Try to find exact or partial matches in the 'romnam' column (case-insensitive)
    country_adm1nm = gdf[gdf['romnam'].str.contains(country_name, case=False, na=False)]

    # If no direct match is found, use fuzzy matching
    if len(country_adm1nm) == 0:
        gdf['distance'] = gdf['romnam'].apply(get_distance)
        matching_rows = gdf[gdf['distance'] < thr]

        if len(matching_rows) > 0:
            print('A direct match was not found. These are the rows closer to your input:')
            print(matching_rows[['romnam', 'distance']])
            print('Please choose a country name from this list or try another country name.')
        else:
            print('No country matching found. Please try a different country name.')

        AOI = None
    else:
        # Dissolve the matched rows into a single geometry
        AOI = country_adm1nm.dissolve()
        print(f'The country "{AOI['romnam'][0]}" was selected from the UN  database.')
        
        
    # Return the result
    if return_gdf:
        return AOI[['romnam', 'iso3cd','maplab','geometry']], gdf
    else:
        return AOI[['romnam', 'iso3cd','maplab','geometry']]


def export_cava_data(ds, years_list, output_folder='./data_input/CAVA_data/', country_name=''):
    """
    Exports selected time slices and rearranged variable data from an xarray dataset.

    Parameters:
    - ds: xarray.Dataset containing the variables and time dimension.
    - years_list: list of years (as integers or strings) to filter the time dimension.
    - output_folder: path to the folder where .npy files will be saved.
    """
    os.makedirs(output_folder, exist_ok=True)

    times_cava = ds['tasmax']['time'].values
    red_times_idx = [i for i, t in enumerate(times_cava) if str(years_list[0]) in str(t)]

    # Save filtered time array
    np.save(os.path.join(output_folder, "time_array.npy"), times_cava[red_times_idx])

    # Export each variable
    for var in ds.keys():
        print(f'Exporting {var}...')
        data = ds[var].values

        # setting proper interval for humidity
        if var == 'hurs':
            if np.max(data) > 1:
                # put data from percentage to the interval [0,1]
                data = data / 100

        rearranged_data = np.moveaxis(data, 0, -1)
        np.save(os.path.join(output_folder, f"{var}_{country_name}.npy"), rearranged_data[:, :, red_times_idx])

    
    # Convert to DataFrame and save
    print("Exporting dataset as DataFrame...")
    df = ds.to_dataframe().reset_index()
    df.to_csv(os.path.join(output_folder, f"cava_dataset_{country_name}.csv"), index=False)


def rasterize_country_from_geodataframe(ds, AOI, country_name, output_dir='./data_input/'):
    """
    Rasterizes a country geometry using spatial metadata from an xarray dataset and a GeoDataFrame.

    Parameters:
    - ds: xarray.Dataset containing 'tasmax', 'latitude', 'longitude'
    - AOI: geopandas.GeoDataFrame containing the country geometry
    - country_name: string, name of the country for file naming
    - output_dir: directory to save the output raster
    """
    os.makedirs(output_dir, exist_ok=True)

    # Extract geometry from GeoDataFrame
    AOI_geom = AOI.geometry

    # Extract spatial metadata from dataset
    try:
        latitudes = ds['pr']['latitude'].values
        longitudes = ds['pr']['longitude'].values
    except:
        latitudes = ds['pr']['lat'].values
        longitudes = ds['pr']['lon'].values
    data_shape = ds['pr'].shape[1:]  # Assuming [time, lat, lon]

    # Define transform and CRS
    transform = from_bounds(min(longitudes), min(latitudes), max(longitudes), max(latitudes),
                            data_shape[1], data_shape[0])
    crs = CRS.from_epsg(4326)

    # Rasterize the country geometry
    country_raster = rasterize(
        [(geom, 1) for geom in AOI_geom],
        out_shape=data_shape,
        transform=transform,
        fill=0,
        dtype=np.uint8
    )

    # Save the rasterized country
    output_path = os.path.join(output_dir, f'{country_name}_rasterized.tif')
    meta = {
        'driver': 'GTiff',
        'height': data_shape[0],
        'width': data_shape[1],
        'count': 1,
        'dtype': np.uint8,
        'crs': crs,
        'transform': transform
    }

    with rasterio.open(output_path, 'w', **meta) as dest:
        dest.write(country_raster, 1)

    print(f"Rasterized country saved to {output_path}")


def clip_raster_with_geometry(raster_path, AOI, country_name, output_dir='./data_input/'):
    """
    Clips a raster using the geometry from a GeoDataFrame and saves the result.

    Parameters:
    - raster_path: path to the input raster file (e.g., DEM)
    - AOI: geopandas.GeoDataFrame containing the geometry to clip with
    - country_name: name of the country for naming the output file
    - output_dir: directory to save the clipped raster
    """
    os.makedirs(output_dir, exist_ok=True)
    AOI_geom = AOI.geometry

    with rasterio.open(raster_path) as src:
        clipped_data, clipped_transform = mask(src, AOI_geom, crop=True)
        clipped_meta = src.meta.copy()
        clipped_meta.update({
            "driver": "GTiff",
            "height": clipped_data.shape[1],
            "width": clipped_data.shape[2],
            "transform": clipped_transform
        })

    output_path = os.path.join(output_dir, f"{country_name}_elevation_clipped.tif")
    with rasterio.open(output_path, 'w', **clipped_meta) as dst:
        dst.write(clipped_data)

    print(f"Clipped raster saved to {output_path}")
    return output_path

def process_raster(INPUT_raster_path, ref_raster_path, output_clipped_path, ref_geom, plot=False, 
                   title = [], colorbar_label = [], cmap = 'viridis', nodata = 0):
    """
    Process a raster by clipping, resampling, and optionally plotting it.
    
    Parameters:
    - INPUT_raster_path: Path to the input raster file.
    - ref_raster_path: Path to the reference raster (e.g., temperature raster).
    - output_clipped_path: Path to save the final clipped raster.
    - ref_geom: Geometry to clip the raster.
    - plot: Boolean flag to plot the result (default: False).
    """
    
    # Load the reference raster (e.g., temperature raster) to get its resolution, extent, and CRS
    with rasterio.open(ref_raster_path) as ref_raster:
        ref_transform = ref_raster.transform
        ref_crs = ref_raster.crs
        ref_shape = ref_raster.shape  # Get the shape (height, width) of the reference raster
        ref_res = ref_raster.res      # Get the pixel resolution of the reference raster

    # Step 1: Load the DEM raster and clip it using the country geometry (bounding box)
    with rasterio.open(INPUT_raster_path) as INPUT_raster:
        raster_clipped, raster_transform = mask(INPUT_raster, ref_geom.envelope, crop=True)
        raster_nodata = INPUT_raster.nodata  # Get the NoData value from the original DEM
        raster_clipped = np.where(raster_clipped == raster_nodata, 0, raster_clipped)  # Replace NoData with 0

        # Step 2: Resample the DEM to match the resolution and extent of the reference raster
        raster_resampled = np.empty(ref_shape, raster_clipped.dtype)  # Create an empty array for resampling
        reproject(
            source=raster_clipped,
            destination=raster_resampled,
            src_transform=raster_transform,
            src_crs=INPUT_raster.crs,
            dst_transform=ref_transform,
            dst_crs=ref_crs,
            resampling=Resampling.med  # Bilinear interpolation for resampling
        )

        # Step 3: Update metadata to match the reference raster's characteristics
        ratser_meta = INPUT_raster.meta.copy()
        ratser_meta.update({
            "driver": "GTiff",
            "height": ref_shape[0],
            "width": ref_shape[1],
            "transform": ref_transform,
            'dtype': np.float32,  # Output data type (unsigned 8-bit)
            "crs": ref_crs,
            "nodata": nodata  # Specify the NoData value
        })

        # Step 4: Save the resampled DEM to a new file
        with rasterio.open(output_clipped_path, 'w', **ratser_meta) as dst:
            dst.write(raster_resampled, 1)  # Write the resampled DEM to disk

    # Step 5: Load the resampled DEM and clip it using the country geometry
    with rasterio.open(output_clipped_path) as resampled_raster:
        resampled_clipped, _ = mask(resampled_raster, ref_geom, crop=True, nodata = -9999)
        ratser_meta = resampled_raster.meta.copy()
        resampled_clipped = resampled_clipped[0]  # Extract the first band (DEM data)

    # Step 6: Save the clipped DEM to a new file
    with rasterio.open(output_clipped_path, 'w', **ratser_meta) as dst:
        dst.write(resampled_clipped, 1)  # Write the clipped DEM

    # Step 7: Optionally, mask -9999 values and plot the final result
    masked_raster = np.ma.masked_equal(resampled_clipped, -9999)

    print('Unique values:', np.unique(masked_raster)[:10])
    
    if plot:
        vmax, vmin = np.max(masked_raster), np.min(masked_raster)
        plt.imshow(masked_raster, vmax=vmax, vmin=vmin, cmap=cmap)  # Use terrain colormap
        plt.colorbar(label=colorbar_label)
        plt.title(title)
        plt.show()

    return masked_raster