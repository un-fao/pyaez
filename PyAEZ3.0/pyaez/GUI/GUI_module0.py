'''import supporting libraries'''
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
try:
    from osgeo import gdal
except:
    import gdal
import sys
import geopandas as gpd
import cavapy
import xarray as xr
import rasterio
from rasterio.transform import from_bounds
from rasterio.crs import CRS
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.features import rasterize
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling
import streamlit as st
import contextlib
import io

os.environ['PROJ_LIB'] = '/opt/conda/envs/pyAEZ/share/proj'
gdal.UseExceptions()

# Add the module0 directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), r'/workspaces/Python-AEZ-FAO')))

# Now you can import your module
from module0 import pre_proc_utilities
import pyaez


# Title
st.title("🌍 Module 0 - Data Retrieval and Pre-processing")

# Section: Climate Model Configuration
st.subheader("🧪 Climate Model Parameters")

# RCP selection
rcp = st.selectbox("🌡️ Select RCP", options=["rcp26", "rcp85"])

# GCM selection
gcm = st.selectbox("🌀 Select GCM", options=["MOHC", "NCC", "MPI"])

# RCM selection
rcm = st.selectbox("🌍 Select RCM", options=["ICTP", "REMO"])

# Year selection
st.subheader("📅 Select Years")

year_mode = st.radio("Choose year input mode:", ["Single Year", "Year Range"])

if year_mode == "Single Year":
    single_year = st.number_input("Enter a year", min_value=1900, max_value=2100, step=1)
    years_list = [single_year]
else:
    start_year = st.number_input("Start Year", min_value=1900, max_value=2100, step=1)
    end_year = st.number_input("End Year", min_value=start_year, max_value=2100, step=1)
    years_list = list(range(start_year, end_year + 1))

# Historical data toggle
historical = st.checkbox("📜 Include Historical Data", value=True)

# Display selected values in a styled box
st.markdown("### 🔎 Selected Parameters")
st.info(f"""
- **RCP**: {rcp}  
- **GCM**: {gcm}  
- **RCM**: {rcm}  
- **Historical Data**: {'Yes' if historical else 'No'}  
- **Years Selected**: {', '.join(map(str, years_list))}
""")


# Section: Country Selection
st.subheader("🗺️ Country Selection")

# Path to the UN shapefile
UN_maps = r'/workspaces/Python-AEZ-FAO/data_input/UN_world_shapefile/BNDA_A1.shp'

# Read the shapefile into a GeoDataFrame
gdf = gpd.read_file(UN_maps)

# Extract country names from the 'romnam' column
country_list = sorted(gdf['romnam'].dropna().unique())

# Add a selection box for country
country_name = st.selectbox("🌐 Select a Country", options=country_list)

# Display selected country
st.success(f"✅ Selected Country: **{country_name}**")

# Try to find exact or partial matches in the 'romnam' column (case-insensitive)
country_adm1nm = gdf[gdf['romnam'].str.contains(country_name, case=False, na=False)]

# Dissolve the matched rows into a single geometry
AOI = country_adm1nm.dissolve()

# Run CAVAPY button
if st.button("🚀 Run CAVAPY"):
    years_up_to = np.max([np.max(years_list)+1, 2007])   
    years_obs = range(1980, np.min([np.max(years_list)+1, 2019]))
    country_name_capitalized = country_name.capitalize()
    bias_correction = np.max(years_list) >= 2007

    # Show the function call parameters
    st.markdown("### 🧾 CAVAPY Function Call")
    st.code(f"""cavapy.get_climate_data(
    country="{country_name_capitalized}",
    cordex_domain="AFR-22",
    buffer=1,
    rcp="{rcp}",
    gcm="{gcm}",
    rcm="{rcm}",
    num_processes=10,
    years_up_to={years_up_to},
    obs=True,
    bias_correction={bias_correction},
    historical={historical},
    years_obs={list(years_obs)}
)""", language="python")

    # Redirect stdout to capture print statements
    output_buffer = io.StringIO()
    with contextlib.redirect_stdout(output_buffer):
        with st.spinner("⏳ Retrieving climate data..."):
            CAVA_climate_data = cavapy.get_climate_data(
                country=country_name_capitalized, 
                cordex_domain="AFR-22", 
                buffer=1,
                rcp=rcp, 
                gcm=gcm, 
                rcm=rcm, 
                num_processes=10,
                years_up_to=years_up_to, 
                obs=True, 
                bias_correction=bias_correction, 
                historical=historical, 
                years_obs=years_obs
            )

    # Show captured output
    st.markdown("### 📄 CAVAPY Output Log")
    st.code(output_buffer.getvalue())

    st.success("✅ CAVA data fully retrieved!")
    ds = CAVA_climate_data

if st.button("📤 Export CAVA Data"):
    if "ds" in st.session_state:
        with st.spinner("📦 Exporting CAVA data..."):
            pre_proc_utilities.export_cava_data(
                st.session_state.ds, years_list, output_folder='/workspaces/Python-AEZ-FAO/data_input/CAVA_data/',
                country_name=country_name,
            )
        st.success("✅ CAVA data exported successfully!")
    else:
        st.warning("⚠️ Please run CAVAPY first to generate climate data.")

# Button to rasterize country
if st.button("🗺️ Rasterize Country"):
    if "ds" in st.session_state:
        with st.spinner("🔄 Rasterizing country geometry..."):
            pre_proc_utilities.rasterize_country_from_geodataframe(
                ds, AOI, country_name, output_dir='./data_input/'
            )
        st.success("✅ Country rasterization completed!")
    else:
        st.warning("⚠️ Please run CAVAPY first to generate climate data.")

# DEM raster path
DEM_raster_path = r'./data_input/DEM/altmed30s.tif'

# Button to rasterize and process DEM
if st.button("🗻 Rasterize and Process DEM"):
    with st.spinner("🔄 Clipping and resampling DEM..."):
        if "ds" in st.session_state:
            # Clip DEM with country geometry
            pre_proc_utilities.clip_raster_with_geometry(
                DEM_raster_path, AOI, country_name, output_dir='./data_input/'
            )

            # Process and resample DEM
            elv = pre_proc_utilities.process_raster(
                INPUT_raster_path=DEM_raster_path,
                ref_raster_path=f'./data_input/{country_name}_rasterized.tif',
                output_clipped_path=f'./data_input/{country_name}_elevation_resampled_clipped.tif',
                ref_geom=AOI.geometry,
                plot=True,
                title='Resampled DEM',
                colorbar_label='Height (m)',
                cmap='terrain'
            )

            st.success("✅ DEM rasterization and processing completed!")
        else:
            st.warning("⚠️ Please run CAVAPY first to generate climate data.")

# Button to process LULC raster
if st.button("🌾 Process LULC Raster"):
    with st.spinner("🔄 Resampling and clipping LULC raster..."):
        if "ds" in st.session_state:
            lulc = pre_proc_utilities.process_raster(
                INPUT_raster_path='./data_input/LULC/soil_regime_CRUTS32_Hist_8110.tif',
                ref_raster_path=f'./data_input/{country_name}_rasterized.tif',
                output_clipped_path=f'./data_input/{country_name}_lulc_resampled_clipped.tif',
                ref_geom=AOI.geometry,
                plot=True,
                title='LULC',
                colorbar_label='Classes',
                cmap='tab10'
            )
            st.success("✅ LULC raster processed successfully!")
        else:
            st.warning("⚠️ Please run CAVAPY first to generate climate data.")