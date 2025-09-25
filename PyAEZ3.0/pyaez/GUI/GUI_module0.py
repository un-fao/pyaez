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

# Historical data toggle
historical = st.checkbox("📜 Include Historical Data", value=True)

# Display selected values in a styled box
st.markdown("### 🔎 Selected Parameters")
st.info(f"""
- **RCP**: {rcp}  
- **GCM**: {gcm}  
- **RCM**: {rcm}  
- **Historical Data**: {'Yes' if historical else 'No'}
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
