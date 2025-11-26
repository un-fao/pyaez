# PyAEZ v3.0 (2025) — Module 4: Soil Constraints (Part 1)

Extract HWSD v2 Soil Mapping Units (SMUs) and slope classes for a user-defined Area of Interest (AOI) across seven layers (D1–D7), clean attributes, and export layer tables for downstream evaluation.


## Authors / credits

* N. Lakmal Deshapriya. 2016 
* Swun Wunna Htet. 2023 
* Swun Wunna Htet. 2024 (Dec)
* Swun Wunna Htet. 2025 (Apr)
* RutendoTadiwa Mukaratirwa (FAO\_NSLD). 2025 (May)
* Shahla Asgharinia (FAO\_NSLD). 2025 (July)

## Summary

This module clips the HWSD raster to an AOI, extracts unique SMUs, and corresponding soil slope classes, joins them to the HWSD attribute workbook, applies light cleaning and transformations, and writes Topsoil (D1–D3) and Subsoil (D4–D7) Excel workbooks. Optional QC images are saved for quick visual checks.

## Key updates (v3.0)

1. Data source harmonized to the latest HWSD release:https://data.apps.fao.org/catalog/iso/ff5c613c-75bb-46a9-a162-bc728059b465 
2. AOI-based SMU extraction implemented for seven soil layers D1–D7: D1 (0-20 cm), D2 (20-40 cm), D3 (40-60 cm), D4 (60-80 cm), D5 (80-100 cm), D6 (100-150 cm), and D7 (150-200 cm).
3. implemented slope classes which are extracted from the global slope database as https://data.apps.fao.org/catalog//iso/44950f32-a84d-4784-ac50-1148beeb8597

## Inputs

* HWSD2 raster (`HWSD2_RASTER/HWSD2.bil` with `.hdr/.prj/.stx`)
* AOI boundary (shape file)
* HWSD attribute workbook (`HWSD2_LAYERS.xlsx`, sheet: `HWSD2_LAYERS_2`)
* SOILFER.SLOPE-MED.tif 

## Outputs

* `Top_Soil_Layers.xlsx` (D1)
* `Sub_Soil_Layers.xlsx` (D2–D7)
* `Raster_and_AOI_Boundary_Check.png` (optional QC)
* `Clipped_Raster_Check.png` (optional QC)
* "Clipped_Slope_Check.png" (optional QC)

## Dependencies

* from pathlib import Path (reading directory)
* import os, sys (reading directory)
* rasterio (raster IO, masking)
* geopandas (vector IO, CRS handling)
* pandas, numpy (tabular processing)
* shapely (geometry mapping for masks)
* matplotlib (QC plots; Agg backend)
* openpyxl (Excel writer)
* from rasterio.mask import mask
* from shapely.geometry import mapping
* from rasterio.warp import reproject, Resampling
* import warnings
* from matplotlib import colors

## Install (tip)

```bash
# geospatial stack (conda-forge)
conda install -c conda-forge geopandas rasterio shapely gdal proj
# general python stack (pip)
pip install pandas numpy matplotlib openpyxl
```

## Configuration

Set paths under “User inputs” in the script.

```python
try:
    WORKING_DIR = Path(__file__).resolve().parent
except NameError:
    WORKING_DIR = Path.cwd()

# -----------------------------
# User inputs
# -----------------------------
RASTER_PATH = WORKING_DIR / "hwsd2_raster" / "HWSD2.bil"
EXCEL_PATH  = WORKING_DIR / "HWSD2_LAYERS_2.xlsx"   # new HWSD layer file
SHAPE_PATH  = WORKING_DIR / "Tanzania" / "Tanzania.shp"
SLOPE_TIF   = WORKING_DIR / "SOILFER.SLOPE-MED.tif"  # categorical slope classes 0..10

# -----------------------------
# Output files
# -----------------------------
OUT_TOPSOIL       = WORKING_DIR / "Top_Soil_Layers.xlsx"
OUT_SUBSOIL       = WORKING_DIR / "Sub_Soil_Layers.xlsx"
QC_RASTER_AOI_PNG = WORKING_DIR / "Raster_and_AOI_Boundary_Check.png"
QC_CLIP_PNG       = WORKING_DIR / "Clipped_Raster_Check.png"
QC_SLOPE_PNG      = WORKING_DIR / "Clipped_Slope_Check.png"
OUT_SLP_CSV       = WORKING_DIR / "SMU_SLOPE_CLASSES.csv"

```

## Workflow

### 1) Verify inputs

The script checks that all required inputs exist— the global 30-arc-second (~1 km) HWSD v2 raster of Soil Mapping Unit (SMU) IDs, the AOI shapefile (your area of interest), global slope map and the HWSD soil layers Excel file— and raises a clear error if any are missing.

### 2) Load raster and AOI; align CRS

* Read HWSD raster to get CRS and nodata (fallback `-9999` if absent).
* Load AOI with GeoPandas and reproject to the raster CRS.
* Save a downsampled raster + AOI overlay for a visual check.

> CRS caution — `OGC:CRS84` vs `EPSG:4326`
> `OGC:CRS84` uses lon,lat axis order; `EPSG:4326` is WGS-84 but officially lat,lon. Ensure all layers share the same CRS before masking.

### 3) Clip raster to AOI

Use `rasterio.mask.mask(..., crop=True, nodata=...)` to produce a clipped raster and transform. Optionally save the clipped raster preview.

### 4) Extract unique smus

Mask nodata and compute `np.unique` on the clipped array to get SMUs present in the AOI.

### 5) Load and filter attributes

Read `HWSD2_LAYERS.xlsx` (sheet `HWSD2_LAYERS`) and filter rows where `HWSD2_SMU_ID` is in the AOI’s SMU set.

### 6) Clean and subset

* `fillna(0)` for numeric gaps
* remove `TEXTURE_USDA == 0`
* for each `(HWSD2_SMU_ID, LAYER)`, keep the row with max `SHARE` (dominant component)

### 7) Texture labels
Map each numeric `TEXTURE_USDA` code (1–13) to its USDA texture class (the script replaces the numeric code with the class name during Step 7).

| Code | USDA texture class |
| ---- | ------------------ |
| 1    | Clay (heavy)       |
| 2    | Silty clay         |
| 3    | Clay (light)       |
| 4    | Silty clay loam    |
| 5    | Clay loam          |
| 6    | Silt               | 
| 7    | Silt loam          |
| 8    | Sandy clay         | 
| 9    | Loam               | 
| 10   | Sandy clay loam    | 
| 11   | Sandy loam         |
| 12   | Loamy sand         | 
| 13   | Sand               | 


### 8) Transform flags

Convert `ADD_PROP` to binary where value `3` becomes `1` and others `0`.

### 9) Final schema

Reorder and rename to the compact schema used downstream.

```
CODE, TXT, OC, pH, TEB, BS, CEC_soil, CEC_clay, RSD,
SPR, SPH, ROOTS, IL, DRG, ESP, EC, CCB, GYP, GRC, VSP, LAYER, SHARE
```

### 10) Export Excel workbooks

Write two multi-sheet files:

* `Top_Soil_Layers.xlsx` with D1
* `Sub_Soil_Layers.xlsx` with D2–D7

Each sheet contains the standardized columns for its layer.


## Troubleshooting

* Empty clip or very few SMUs: check CRS alignment and AOI validity. If needed, clean minor topology issues with `buffer(0)` and `dissolve()`.
* NoData handling: if the source raster lacks nodata, the script uses `-9999`. Confirm downstream tools honor it.
* Texture mapping: rows with `TEXTURE_USDA == 0` are dropped by design. Adjust Step 6 if you must retain them.




