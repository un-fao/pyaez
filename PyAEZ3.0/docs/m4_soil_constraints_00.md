# PyAEZ v3.0 (2025) — Module 4: Soil Constraints (Part 1)

Extract HWSD v2 Soil Mapping Units (SMUs) for a user-defined Area of Interest (AOI) across seven soil layers (D1–D7), clean attributes, and export ready-to-use layer tables for downstream Soil Constraints evaluation.

## Authors / Credits

2016: N. Lakmal Deshapriya

2023: Swun Wunna Htet

2024 (Dec): Swun Wunna Htet

2025 (Apr): Swun Wunna Htet

2025 (May): RutendoTadiwa Mukaratirwa (FAO_NSLD)

2025 (July): Shahla Asgharinia (FAO_NSLD)

## Summary

This module clips the HWSD raster to your AOI, extracts unique SMUs, joins them to the HWSD attribute workbook, applies light cleaning/transformations, and writes Topsoil (D1–D3) and Subsoil (D4–D7) Excel workbooks. Optional QC images are saved for fast visual checks.

Key updates (v3.0)

Data source harmonized to the latest HWSD release:
https://data.apps.fao.org/catalog/iso/ff5c613c-75bb-46a9-a162-bc728059b465

AOI-based SMU extraction implemented for seven layers (D1–D7).

## Inputs

HWSD2 raster (HWSD2_RASTER/HWSD2.bil + sidecars .hdr/.prj/.stx)

AOI boundary (vector file with a defined CRS; e.g., Shapefile/GPKG/GeoJSON)

HWSD attribute workbook (HWSD2_LAYERS.xlsx, sheet: HWSD2_LAYERS)

## Outputs

Top_Soil_Layers.xlsx (D1–D3)

Sub_Soil_Layers.xlsx (D4–D7)

Raster_and_AOI_Boundary_Check.png (optional QC)

Clipped_Raster_Check.png (optional QC)

## Dependencies

rasterio (raster IO, masking)

geopandas (vector IO, CRS handling)

pandas, numpy (tabular processing)

shapely (geometry mapping for masks)

matplotlib (QC plots; Agg backend for headless)

openpyxl (Excel writer)

Tip (Windows/macOS): install geospatial stack via conda-forge, then pip for the rest. 

### Install (tip)                                          
```bash
conda install -c conda-forge geopandas rasterio shapely gdal proj
pip install pandas numpy matplotlib openpyxl```

# PyAEZ v3.0 — Module 4: Soil Constraints (Part 1)      <!-- H1: one per page -->

## Authors / Credits                                       <!-- H2 -->
- 2016: N. Lakmal Deshapriya
- 2023: Swun Wunna Htet
- 2025 (July): Shahla Asgharinia (FAO_NSLD)

## Summary                                                 <!-- H2 -->
Extract HWSD v2 SMUs for an AOI across D1–D7...

## Key updates (v3.0)                                      <!-- H2 -->
1. Data source harmonized…
2. AOI-based SMU extraction…

## Inputs                                                  <!-- H2 -->
- HWSD2 raster (.bil + sidecars)
- AOI boundary (vector; CRS defined)
- HWSD attributes Excel

## Outputs                                                 <!-- H2 -->
- Top_Soil_Layers.xlsx (D1–D3)
- Sub_Soil_Layers.xlsx (D4–D7)

## Dependencies                                            <!-- H2 -->
`rasterio`, `geopandas`, `pandas`, `numpy`, `shapely`, `matplotlib`, `openpyxl`

### Install (tip)                                          <!-- H3 -->
```bash
conda install -c conda-forge geopandas rasterio shapely gdal proj
pip install pandas numpy matplotlib openpyxl```

