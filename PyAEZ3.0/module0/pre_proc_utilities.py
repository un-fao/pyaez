"""
PyAEZ version 3.0 (April 2026)
This function library provides support for the execution of Module 0
related to the preparation of input data for Module 1.
6 June 2025: Dario Spiller, initialization of beta version

Modifications:
"""

import pandas as pd
import geopandas as gpd
import os
import numpy as np
import rasterio
import pickle
import matplotlib.pyplot as plt
from rasterio.transform import from_bounds
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling
from typing import Dict
from Levenshtein import distance as lev_distance  
import ipywidgets as widgets
import tkinter as tk
from tkinter import ttk
from pathlib import Path
import sys
from datetime import datetime
import warnings
warnings.filterwarnings("always")

try:
    from osgeo import gdal
except:
    import gdal

from shapely.errors import ShapelyDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="Geometry is in a geographic CRS.*"
    )
    
gdal.UseExceptions()

def set_working_folder():
    """
    Automatically sets the working directory to the project root.

    This function assumes the current Jupyter notebook is located inside the
    `tutorials/` subfolder of the project directory:

        project_root/
            tutorials/
                notebook.ipynb

    The function:
    - Detects the folder of the currently running notebook (using Path().resolve())
    - Moves one level up to reach the project root (`Python-AEZ-FAO`)
    - Changes the working directory to the project root
    - Appends the project root to sys.path for module imports

    Returns
    -------
    Path
        The resolved path of the project root directory.
    """

    # Directory where the notebook is located
    nb_dir = Path().resolve()

    # The project root is the parent folder of tutorials/
    project_root = nb_dir.parent

    # Set working directory
    os.chdir(project_root)
    sys.path.append(str(project_root))

    print(f"Notebook directory : {nb_dir}")
    print(f"Project root set to: {project_root}")

    return None
    
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
        print(f'The country "{AOI["romnam"][0]}" was selected from the UN database.')
        
        
    # Return the result
    if return_gdf:
        return AOI[['romnam', 'iso3cd','maplab','geometry']], gdf
    else:
        return AOI[['romnam', 'iso3cd','maplab','geometry']]


'''    
def choose_country(thr=3, return_gdf=False):
    UN_maps = r'data_input/UN_world_shapefile/BNDA_A1.shp'
    gdf = gpd.read_file(UN_maps)

    country_list = sorted(gdf['romnam'].dropna().unique())

    # ---- Tkinter modal dialog ----
    root = tk.Tk()
    root.title("Select Country")
    root.geometry("300x150")

    ttk.Label(root, text="Choose a country:").pack(pady=10)

    var = tk.StringVar()
    dropdown = ttk.Combobox(root, textvariable=var, values=country_list, state="readonly")
    dropdown.pack(pady=5)

    # storage for selected AOI
    result = {"AOI": None}

    def on_select():
        selected = var.get()
        AOI = get_UN_country(selected, thr=thr, return_gdf=return_gdf)
        result["AOI"] = AOI
        root.destroy()

    ttk.Button(root, text="Select", command=on_select).pack(pady=10)

    # Make window modal (blocks until destroyed)
    root.mainloop()

    return result["AOI"], result["AOI"]['romnam'][0]
'''
def choose_country(thr=3, return_gdf=False):
    import tkinter as tk
    from tkinter import ttk, messagebox
    import geopandas as gpd

    UN_maps = r'data_input/UN_world_shapefile/BNDA_A1.shp'
    gdf = gpd.read_file(UN_maps)

    country_list = sorted(gdf['romnam'].dropna().unique())

    # ---- Tkinter modal dialog ----
    root = tk.Tk()
    root.title("Select Country")
    root.geometry("300x220")
    root.resizable(False, False)

    ttk.Label(root, text="Choose a country:").pack(pady=(10, 2))

    var_country = tk.StringVar()
    dropdown = ttk.Combobox(
        root, textvariable=var_country, values=country_list, state="readonly"
    )
    dropdown.pack(pady=5)

    # ---- Buffer input ----
    ttk.Label(root, text="Buffer (degrees, positive):").pack(pady=(10, 2))

    var_buffer = tk.StringVar(value="0.0")
    buffer_entry = ttk.Entry(root, textvariable=var_buffer)
    buffer_entry.pack(pady=5)

    # storage for selected AOI
    result = {"AOI": None, "buffer": None}

    def on_select():
        try:
            buffer_deg = float(var_buffer.get())
            if buffer_deg < 0:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid input",
                "Buffer must be a positive number (degrees)."
            )
            return

        selected = var_country.get()
        if not selected:
            messagebox.showerror(
                "Selection missing",
                "Please select a country."
            )
            return

        AOI = get_UN_country(
            selected,
            thr=thr,
            return_gdf=return_gdf
        )

        result["AOI"] = AOI
        result["buffer"] = buffer_deg
        root.destroy()

    ttk.Button(root, text="Select", command=on_select).pack(pady=15)

    # Make window modal
    root.mainloop()

    return (
        result["AOI"],
        result["AOI"]["romnam"].iloc[0],
        result["buffer"]
    )    


def choose_scenario():
    """
    Opens a Tkinter dialog where the user selects:
    year, rcp, gcm, rcm, historical, management.
    Returns a dictionary with all selected values.
    """

    root = tk.Tk()
    root.title("Select Scenario Parameters")
    root.geometry("350x400")

    # -------------------
    # Dropdown definitions
    # -------------------

    # YEAR
    ttk.Label(root, text="Select Year:").pack(pady=5)
    year_var = tk.StringVar()
    year_list = [str(y) for y in range(1980, 2100 + 1)]
    ttk.Combobox(root, textvariable=year_var, values=year_list, state="readonly").pack()

    # RCP
    ttk.Label(root, text="Select RCP:").pack(pady=5)
    rcp_var = tk.StringVar()
    rcp_list = ["rcp26", "rcp85"]
    ttk.Combobox(root, textvariable=rcp_var, values=rcp_list, state="readonly").pack()

    # GCM
    ttk.Label(root, text="Select GCM:").pack(pady=5)
    gcm_var = tk.StringVar()
    gcm_list = ["MOHC", "NCC", "MPI"]
    ttk.Combobox(root, textvariable=gcm_var, values=gcm_list, state="readonly").pack()

    # RCM
    ttk.Label(root, text="Select RCM:").pack(pady=5)
    rcm_var = tk.StringVar()
    rcm_list = ["Reg", "REMO"]
    ttk.Combobox(root, textvariable=rcm_var, values=rcm_list, state="readonly").pack()

    # HISTORICAL
    ttk.Label(root, text="Historical data (1980-2005):").pack(pady=5)
    historical_var = tk.StringVar()
    hist_list = ["True", "False"]
    ttk.Combobox(root, textvariable=historical_var, values=hist_list, state="readonly").pack()

    # MANAGEMENT
    ttk.Label(root, text="Management:").pack(pady=5)
    management_var = tk.StringVar()
    management_list = ["HI", "LI"]
    ttk.Combobox(root, textvariable=management_var, values=management_list, state="readonly").pack()

    # Storage for results
    result = {"values": None}

    # -------------------
    # Confirm button
    # -------------------

    def on_confirm():
        result["values"] = {
            "year": int(year_var.get()),
            "rcp": rcp_var.get(),
            "gcm": gcm_var.get(),
            "rcm": rcm_var.get(),
            "historical": (historical_var.get() == "True"),
            "management": management_var.get()
        }
        root.destroy()

    ttk.Button(root, text="Confirm Selection", command=on_confirm).pack(pady=20)

    # Block until user closes dialog
    root.mainloop()

    print('The chosen values for the current simulation are:', result["values"])
    
    return result["values"]




def export_cava_data(ds, year, output_folder='./data_input/CAVA_data/', country_name=''):
    """
    Exports selected time slices and rearranged variable data from an xarray dataset.

    Parameters:
    - ds: xarray.Dataset containing the variables and time dimension.
    - year: year (as integer or string) to filter the time dimension.
    - output_folder: base folder where output subfolder will be created.
    - country_name: optional name used in the subfolder.
    """

    # --- Create timestamp ---
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # numerical up to seconds

    # --- Build subfolder name ---
    if country_name:
        subfolder_name = f"{country_name}_{timestamp}"
    else:
        subfolder_name = timestamp

    # --- Final output path ---
    final_output = os.path.join(output_folder, subfolder_name)
    os.makedirs(final_output, exist_ok=True)

    # --- Select time indices for given year ---
    times_cava = ds['tasmax']['time'].values
    red_times_idx = [i for i, t in enumerate(times_cava) if str(year) in str(t)]

    # --- Save filtered time array ---
    np.save(os.path.join(final_output, "time_array.npy"), times_cava[red_times_idx])

    # --- Export each variable ---
    for var in ds.keys():
        print(f'Exporting {var}...')
        data = ds[var].values

        # Normalize humidity to [0,1] if needed
        if var == 'hurs':
            if np.max(data) > 1:
                data = data / 100

        # Rearrange dimensions and save
        rearranged_data = np.moveaxis(data, 0, -1)
        np.save(os.path.join(final_output, f"{var}.npy"),
                rearranged_data[:, :, red_times_idx])

    print(f"\n✅ Export complete. Files saved in: {final_output}")
    return final_output

def rasterize_country_from_geodataframe(ds, AOI, country_name, output_dir='./data_input/'):
    """
    Rasterizes a country geometry using spatial metadata from an xarray dataset and a GeoDataFrame.
    Also displays a simple plot to verify the rasterization.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Extract geometry
    AOI_geom = AOI.geometry

    # Extract spatial metadata
    try:
        latitudes = ds['pr']['latitude'].values
        longitudes = ds['pr']['longitude'].values
    except:
        latitudes = ds['pr']['lat'].values
        longitudes = ds['pr']['lon'].values

    data_shape = ds['pr'].shape[1:]   # (lat, lon)

    # Define transform and CRS
    transform = from_bounds(min(longitudes), min(latitudes),
                            max(longitudes), max(latitudes),
                            data_shape[1], data_shape[0])
    crs = CRS.from_epsg(4326)

    # Rasterize geometry (1 = inside AOI, 0 = outside)
    country_raster = rasterize(
        [(geom, 1) for geom in AOI_geom],
        out_shape=data_shape,
        transform=transform,
        fill=0,
        dtype=np.uint8
    )

    # Save raster
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

    # ✅ Simple verification plot
    plt.figure(figsize=(6, 6))
    plt.imshow(country_raster, cmap='Greens')
    plt.title(f"Rasterized Country Mask: {country_name}")
    plt.xlabel("Longitude pixels")
    plt.ylabel("Latitude pixels")
    plt.grid(False)
    plt.show()

    return country_raster



def clip_raster_with_geometry(raster_path, AOI, country_name, buffer, output_dir='./data_input/'):
    """
    Clips a raster using the geometry from a GeoDataFrame, saves the result, 
    and shows a plot of the clipped raster.
    """

    os.makedirs(output_dir, exist_ok=True)
    AOI_geom = AOI.geometry.envelope.buffer(buffer)

    # --- Clip raster ---
    with rasterio.open(raster_path) as src:
        clipped_data, clipped_transform = mask(src, AOI_geom, crop=True, nodata = -9999)
        clipped_meta = src.meta.copy()
        clipped_meta.update({
            "driver": "GTiff",
            "height": clipped_data.shape[1],
            "width": clipped_data.shape[2],
            "transform": clipped_transform
        })

    # --- Save clipped raster ---
    output_path = os.path.join(output_dir, f"{country_name}_elevation_clipped.tif")
    with rasterio.open(output_path, 'w', **clipped_meta) as dst:
        dst.write(clipped_data)

    print(f"Clipped raster saved to {output_path}")
    
    # --- Plot clipped raster ---
    masked_raster = np.ma.masked_equal(clipped_data[0], -9999)
    
    plt.figure(figsize=(7, 6))

    # If raster has multiple bands, show band 1
    plt.imshow(masked_raster, cmap='terrain')
    plt.title(f"Clipped Elevation Raster - {country_name}")
    plt.xlabel("Column index")
    plt.ylabel("Row index")
    plt.colorbar(label="Elevation")
    plt.grid(False)
    plt.show()

    return output_path
    
def process_raster(INPUT_raster_path, ref_raster_path, output_clipped_path, ref_geom, buffer, plot=False, 
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
        raster_clipped, raster_transform = mask(INPUT_raster, ref_geom.envelope.buffer(buffer), crop=True)
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
            resampling=Resampling.average  # Bilinear interpolation for resampling
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
        resampled_clipped, _ = mask(resampled_raster, ref_geom.envelope.buffer(buffer), crop=True, nodata = -9999)
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

def choose_crop(crop_param_filename):
    """
    Opens a dropdown list of available crops extracted from the Excel parameter file.
    Returns the selected crop name as a string.
    """

    # Load crop table
    df = pd.read_excel(crop_param_filename, sheet_name='parameters')

    # Extract unique crop names
    available_crops = df[df['Crop_name'].notna()]['Crop_name'].unique()

    # ------- Tkinter GUI -------
    root = tk.Tk()
    root.title("Select Crop")
    root.geometry("350x150")

    ttk.Label(root, text="Choose a crop:").pack(pady=10)

    crop_var = tk.StringVar()
    cb = ttk.Combobox(
        root, 
        textvariable=crop_var, 
        values=list(available_crops),
        state="readonly",
        width=40
    )
    cb.pack(pady=5)

    # Storage for result
    result = {"crop": None}

    def on_select():
        result["crop"] = crop_var.get()
        root.destroy()

    ttk.Button(root, text="Select", command=on_select).pack(pady=15)

    # Block until user makes a selection
    root.mainloop()

    print(f'These are the parameters selected for {result["crop"]}:\n')
    display(df[df['Crop_name'] == result["crop"]])
    
    return result["crop"]

def select_from_dropdown(options, title="Select an option", label="Choose one:"):
    """
    Display a modal Tkinter dropdown list and return the selected value.
    `options` must be a list of tuples: (value_to_return, text_to_display)
    """
    root = tk.Tk()
    root.title(title)
    root.geometry("450x180")

    ttk.Label(root, text=label).pack(pady=10)

    var = tk.StringVar()

    # Dropdown items only show the text part
    dropdown = ttk.Combobox(root,
                            textvariable=var,
                            values=[text for _, text in options],
                            state="readonly",
                            width=50)
    dropdown.pack(pady=5)

    result = {"value": None}

    def on_select():
        idx = dropdown.current()
        result["value"] = options[idx][0]   # return the real underlying value
        root.destroy()

    ttk.Button(root, text="Select", command=on_select).pack(pady=10)
    root.mainloop()

    return result["value"]

def get_agro_climatic_constraints_array_full_appendices_M3(
    excel_path: str,
    crop_name: str,
    input_level: str,
    condition_type: str = 'rainfed',
    growing_length_period: int = 0,
    crop_column: str = "Common name",
    glp_column: str = "Growth Cycle",
    drm_column: str = "Dormancy"
) -> dict:  
    """
    Extract agro‑climatic constraint arrays from the official Appendix 5-1 tables.

    This function:
    - Selects the correct three sheets depending on `condition_type`
      * rainfed    → A5‑1.1 (gte20), A5‑1.2 (lt10), A5‑1.5 (lgpt10)
      * irrigated  → A5‑1.3 (gte20), A5‑1.4 (lt10), A5‑1.5 (lgpt10)
    - Cleans and harmonizes each sheet:
      * drops extra structural columns
      * renames crop‑specific metadata columns (Dormancy, Growth Cycle)
      * forward-fills metadata within each crop group
      * filters by input level (except for sheet 1.5)
      * filters by crop name (partial case‑insensitive match)
      * filters by growing length period if specified
    - Handles multiple matching rows by asking the user to choose one.
      For A5‑1.5 it returns 1 row; for others it returns 4 rows starting
      from the selected index.
    - Extracts only the agro‑climatic constraint columns, harmonizes
      their names, and returns them as DataFrames.

    Returns
    -------
    dict
        A dictionary with three keys:
        - "gte20": DataFrame of frost‑free constraints for ≥20% probability
        - "lt10":  DataFrame for <10% probability
        - "lgpt10": DataFrame for >10 months frost‑free period
    """

    # ----------
    # Select sheets based on condition type
    # ----------
    if condition_type == "rainfed":
        sheet_map = {
            "gte20": ["A5-1.1", "Agro-climatic constraints for rain-fed conditions when mean annual temperature > 20 degree"],
            "lt10": ["A5-1.2", "Agro-climatic constraints for rain-fed conditions when mean annual temperature < 10 degree"],
            "lgpt10": ["A5-1.5", "Agro-climatic constraints of early and late frost"]
        }
    elif condition_type == "irrigated":
        sheet_map = {
            "gte20": ["A5-1.3", "Agro-climatic constraints for irrigated conditions when mean annual temperature > 20 degree"],
            "lt10": ["A5-1.4", "Agro-climatic constraints for irrigated conditions when mean annual temperature < 10 degree"],
            "lgpt10": ["A5-1.5", "Agro-climatic constraints of early and late frost"]
        }
    else:
        raise ValueError('condition_type must be either "rainfed" or "irrigated"')

    results: Dict[str, pd.DataFrame] = {}

    # Loop through the three sheets
    for key, (sname, constraint_type) in sheet_map.items():
 
        print('')
        print(constraint_type)
        print('Please select the crop of interest:')
        # ----------
        # Load sheet
        # ----------
        df = pd.read_excel(
            excel_path,
            sheet_name=sname,
            skiprows=3,                       # skip metadata rows
            usecols=lambda c: c not in ['Unnamed: 0']  # ignore Excel index column
        )

        # Drop Excel's third structural column
        df = df.drop(df.columns[2], axis=1)

        # Rename metadata columns by position
        df = df.rename(columns={
            df.columns[1]: drm_column,
            df.columns[2]: glp_column
        })

        # Forward‑fill crop metadata within each crop group
        df[[drm_column, glp_column]] = (
            df.groupby("Common name")[[drm_column, glp_column]].ffill()
        )

        # Remove empty rows
        df = df.dropna(how='all')

        # ----------
        # Filter by input level and constraint type (except sheet 1.5)
        # ----------
        if sname != "A5-1.5":
            level_mask = df["Input level"].str.contains(
                input_level, case=False, na=False
            )
            constraint_type_mask = df["Constraint type"] != 'a'
            df = df[level_mask & constraint_type_mask]

        # Filter by crop name (partial match, case-insensitive)
        name_mask = df[crop_column].str.contains(
            crop_name, case=False, na=False
        )

        filtered_df_name = df[name_mask]
        
        # Filter by LGP if specified
        if growing_length_period != 0:
            glp_mask = filtered_df_name[glp_column] == growing_length_period
            
            # If no rows match, fall back to True mask + warn user
            if not glp_mask.any():
                warnings.warn(
                    "Be aware that no correspondence was found for the requested length "
                    "of the growing period; please select another set of constraints.",
                    UserWarning
                )
                glp_mask = pd.Series(True, index=filtered_df_name.index)
        
        else:
            # Case LGP = 0 → accept all rows
            glp_mask = pd.Series(True, index=filtered_df_name.index)

        filtered_df = filtered_df_name[glp_mask]

        # No matches found
        if filtered_df.empty:
            raise ValueError(
                f"No data found for crop '{crop_name}' with LGP={growing_length_period}"
            )

        # ----------
        # Resolve row selection
        # ----------
        if len(filtered_df) == 1:
            # Unique match → trivial
            selected_row = filtered_df.iloc[[0]]

        elif len(filtered_df) > 1:
            # ---------- Dropdown GUI for disambiguation ----------
            
            print("\n--- Multiple matching rows found ---")
        
            options = []
            for idx, row in filtered_df.iterrows():
                if row["Constraint type"] in ('b', 'e'):
                    text = (f"Row {idx} → {row[crop_column]} | "
                            f"GLP={row[glp_column]} | Input={row['Input level']} | "
                            f"Type={row['Constraint type']}")
                    options.append((idx, text))
        
            selected_idx = select_from_dropdown(
                options,
                title="Select Crop Variant",
                label="Multiple matches found — please select the correct row:"
            )
        
            # A5‑1.5 returns one row; others return 4 rows
            if sname == "A5-1.5":
                selected_row = df.loc[[selected_idx]]
            else:
                selected_row = df.loc[selected_idx:selected_idx + 3]
        

        # ----------
        # Extract constraint columns
        # ----------
        if sname != "A5-1.5":
            # Use column index 4 + columns from index 6 forward
            columns = [df.columns[4]] + list(df.columns[6:])
            new_cols_name = [
                'type','0,29','30,59','60,89','90,119','120,149','150,179',
                '180,209','210,239','240,269','270,299','300,329',
                '330,364','365-','365+'
            ]
        else:
            # A5‑1.5 uses all columns from index 4 forward
            columns = list(df.columns[4:])
            new_cols_name = [
                'type','0,29','30,59','60,89','90,119','120,149','150,179',
                '180,209','210,239','240,269','270,299','300,329',
                '330,364','365,366'
            ]

        #display(selected_row[columns])
        
        # Store final processed array
        results[key] = (
            selected_row[columns]
            .replace("n.a.", 0)
            .replace("NaN", 0)
            .replace(np.nan, 0)
            .rename(columns=dict(zip(columns, new_cols_name)))
        )
        #print('ciao!!')
        #display(results)
        print('These are the contraints you have selected:')
        display(results[key])
    
    for key, df in results.items():
        # Convert all constraint columns except 'type' to numeric
        numeric_cols = df.columns[1:]   # everything except 'type'
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0).astype('int64')

    with pd.HDFStore(f"./data_input/input_module3/M3_constraints_{crop_name}_{input_level}_{condition_type}.h5") as store:
        for key, df in results.items():
            store.put(key, df, format="fixed")
        
    return results


def initialize_clim(
        year: str,
        country_name: str,
        daily: bool = True,
        climate_folder: str = ''
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
    daily : bool, optional
        If True, filter climate arrays for the given year and pass daily data.
        If False, pass full arrays (monthly or climatology).
    
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
    ...     daily=True,
    ... )
    >>> # ClimateRegime object is now ready for NB1 computations.
    """

    mask_value = 0  # pixel value in admin_mask to exclude from the analysis
    country_mask_name = f'{country_name}_rasterized.tif'
    work_dir=os.getcwd()

    ELV_filename = f'{country_name}_elevation_resampled_clipped.tif'
    ELV_raster_path = os.path.join(work_dir, 'data_input', ELV_filename)
    elevation = gdal.Open(ELV_raster_path).ReadAsArray()

    # Validate working directory
    if not os.path.isdir(work_dir):
        raise FileNotFoundError(f"Working directory not found: {work_dir}")

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB1')
    os.makedirs(folder_path, exist_ok=True)

    from pyaez import ClimateRegime

    # Initialize ClimateRegime object
    clim_reg = ClimateRegime.ClimateRegime()

    # Load climate data
    data_input = os.path.normpath(os.path.join(work_dir, climate_folder))
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

    red_times_idx = [i for i,t in enumerate(times) if str(year) in str(t)]
    
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

    clim_reg.country_name = country_name
    clim_reg.year = year
    clim_reg.ref_raster = mask_path
    
    # Path to directory and file
    dir_path = r'./data_input/input_module1'
    file_path = os.path.join(dir_path, 'clim_reg.pkl')
    
    # Create directory if it does not exist
    os.makedirs(dir_path, exist_ok=True)
    
    # Save the pickle file
    with open(file_path, "wb") as f:
        pickle.dump(clim_reg, f)

    return clim_reg

def initialize_aez(
        year: str,
        country_name: str,
        management: str = 'HI',
        daily: bool = True,
        climate_folder: str = ''
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

    country_mask_name = f'{country_name}_rasterized.tif'

    # Set working directory
    work_dir=os.getcwd()
    #work_dir = os.path.dirname(cwd)
    # sys.path.append(work_dir)

    
    # Import your library
    from pyaez import CropSimulation
    from pyaez.UtilitiesCalc import saveRaster, classifyFinalYield

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB2')
    os.makedirs(folder_path, exist_ok=True)

    # Initialize CropSimulation
    aez = CropSimulation.CropSimulation()

    # Load climate data
    data_input = os.path.normpath(os.path.join(work_dir, climate_folder))
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
    
    ELV_filename = f'{country_name}_elevation_resampled_clipped.tif'
    ELV_raster_path = os.path.join(work_dir, 'data_input', ELV_filename)
    elevation = gdal.Open(ELV_raster_path).ReadAsArray()

    # Extract geotransform and bounds
    geotransform = mask_gdal.GetGeoTransform()
    minx = geotransform[0]
    miny = geotransform[3] + geotransform[5] * mask_gdal.RasterYSize
    maxx = geotransform[0] + geotransform[1] * mask_gdal.RasterXSize
    maxy = geotransform[3]

    lat_min, lat_max = miny, maxy

    # Filter times for the given year
    red_times_idx = [i for i, t in enumerate(times) if str(year) in str(t)]
    if len(red_times_idx) == 366:  # Leap year adjustment
        red_times_idx = np.delete(red_times_idx, 59)

    # Configure AEZ objectcrop_param_filename
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
        crop_param_filename = r'./data_input/input_crop_TSUM_parameters_HI.xlsx'
    elif management == 'LI':
        # crop_param_filename = r'./data_input/input_crop_TSUM_parameters_LI.xlsx'
        raise RuntimeError("Low input parameters are yet to be provided")
    else:
        raise NameError("The variable management can be either 'HI' or 'LI'.")

    crop_name = choose_crop(crop_param_filename)

    
    #df = pd.read_excel(crop_param_filename, sheet_name='parameters')
    #
    #
    ## Filter non-NaN crop names
    #available_crops = df[df['Crop_name'].notna()]['Crop_name'].unique()
    #
    #print("Available crops:")
    #for i, crop in enumerate(available_crops, start=1):
    #    print(f"{i}. {crop}")
    #
    #choice = int(input("Select a crop by number: "))
    #crop_name = available_crops[choice - 1]
    #print(f"You selected: {crop_name}")
    
    aez.readCropandCropCycleParameters(file_path = crop_param_filename, 
                                   crop_name = crop_name)

    # Updated Temperature Profile screenign routine
    aez.setupCropSpecificRule(file_path = r'./data_input/crop-specific_rule.xlsx',
                                   crop_name = crop_name.split('_')[0])

    # Load mask and elevation
    mask_path = os.path.join(work_dir, 'data_input', country_mask_name)
    if not os.path.exists(mask_path):
        raise FileNotFoundError(f"Mask file not found: {mask_path}")

    aez.country_name = country_name
    aez.year = year
    aez.ref_raster = mask_path
    
    # Path to directory and file
    dir_path = r'./data_input/input_module2'
    file_path = os.path.join(dir_path, 'aez.pkl')
    
    # Create directory if it does not exist
    os.makedirs(dir_path, exist_ok=True)
    
    # Save the pickle file
    with open(file_path, "wb") as f:
        pickle.dump(aez, f)

    return aez

def initialize_climatic_constraints(
        year: str,
        country_name: str,
        management: str,
        aez,
        climate_folder: str
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
    year : str
        Year selector used to filter the time dimension of climate arrays.
        Example: `"1990"`. If the selected year has 366 daily entries, day index 59
        (Feb 29) is removed to keep 365 days.
    country_name : str
        Country or study-area name (used for file naming and logging).
    management : {'HI','LI'}
        Crop parameter set:
        - 'HI' loads `input_crop_TSUM_parameters_HI.xlsx`,
        - 'LI' currently raises `RuntimeError` (not implemented).
    
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
    ...        ...     management="HI",
    ...     year="1990"
    ... )
    >>> # ClimaticConstraints object is now ready for NB3 computations.
    """

    country_mask_name = f'{country_name}_rasterized.tif'
    
    # Set working directory
    work_dir=os.getcwd()
    
    # Import AEZ library
    from pyaez import ClimaticConstraints

    # Ensure output folder exists
    folder_path = os.path.join(work_dir, 'data_output', 'NB3')
    os.makedirs(folder_path, exist_ok=True)

    # Load mask and elevation
    mask_path = os.path.join(work_dir, 'data_input', country_mask_name)
    mask_gdal = gdal.Open(mask_path)
    mask = mask_gdal.ReadAsArray()
    
    ELV_filename = f'{country_name}_elevation_resampled_clipped.tif'
    ELV_raster_path = os.path.join(work_dir, 'data_input', ELV_filename)
    elevation = gdal.Open(ELV_raster_path).ReadAsArray()

    # Extract geotransform and bounds
    geotransform = mask_gdal.GetGeoTransform()
    lat_min = geotransform[3] + geotransform[5] * mask_gdal.RasterYSize
    lat_max = geotransform[3]

    # Load climate data
    data_input = os.path.normpath(os.path.join(work_dir, climate_folder))
    max_temp = np.load(os.path.join(data_input, 'tasmax.npy'))
    min_temp = np.load(os.path.join(data_input, 'tasmin.npy'))
    precipitation = np.load(os.path.join(data_input, 'pr.npy'))
    rel_humidity = np.load(os.path.join(data_input, 'hurs.npy'))
    wind_speed = np.load(os.path.join(data_input, 'sfcWind.npy'))
    short_rad = np.load(os.path.join(data_input, 'rsds.npy'))
    times = np.load(os.path.join(data_input, 'time_array.npy'))

    
    # Filter times for the given year
    red_times_idx = [i for i, t in enumerate(times) if str(year) in str(t)]
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

    excel_file = r"./data_input/GAEZ4_Appendices.xlsx"

    _ = get_agro_climatic_constraints_array_full_appendices_M3(
        excel_path=excel_file,
        crop_name=aez.crop_name,
        input_level=management,
        condition_type = 'rainfed',
        growing_length_period = aez.cycle_len)
    
    _ = get_agro_climatic_constraints_array_full_appendices_M3(
        excel_path=excel_file,
        crop_name=aez.crop_name,
        input_level=management,
        condition_type = 'irrigated',
        growing_length_period = aez.cycle_len)
    
    if not hasattr(clim_con, 'crop_name'):
        setattr(clim_con, 'crop_name', aez.crop_name)
    else:
        clim_con.crop_name = aez.crop_name
        
    clim_con.input_level = management
    clim_con.country_name = country_name
    clim_con.year = year
        
    with open(rf'./data_input/input_module3/clim_con.pkl', "wb") as f:
        pickle.dump(clim_con, f)
    
    return clim_con