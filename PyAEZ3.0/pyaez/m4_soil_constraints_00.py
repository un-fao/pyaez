"""
PyAEZ v3.0 (2025)
Module 4 — Soil Constraints (Part 1)

Authors / Credits
- 2016: N. Lakmal Deshapriya
- 2023: Swun Wunna Htet
- 2024 (Dec): Swun Wunna Htet
- 2025 (Apr): Swun Wunna Htet
- 2025 (May): RutendoTadiwa Mukaratirwa (FAO_NSLD)
- 2025 (July): Shahla Asgharinia (FAO_NSLD)


Summary
This module extracts Harmonized World Soil Database (HWSD v2) Soil Mapping Units (SMUs)
for a user-defined Area of Interest (AOI) across seven standard soil layers (D1–D7).
It prepares cleaned, layer-wise tables to feed subsequent Soil Constraints evaluation.

Key updates (v3.0)
1) Data source harmonized to the latest HWSD release:
   https://data.apps.fao.org/catalog/iso/ff5c613c-75bb-46a9-a162-bc728059b465
2) Implements AOI-based SMU extraction for seven soil layers (D1–D7).

Inputs
- HWSD2 raster (BIL + sidecars)
- AOI boundary (vector; CRS defined)
- HWSD attribute workbook (Excel)

Outputs
- Cleaned per-layer tables for D1–D7 (e.g., Topsoil D1–D3, Subsoil D4–D7)
- Optional QC plots for raster/AOI overlay and clipped raster
"""

from pathlib import Path
import os
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import mapping
from openpyxl import Workbook

# Dependencies (why each is used)
# - pathlib.Path: safe, cross-platform file paths.
# - rasterio (+ rasterio.mask): read/clip HWSD .bil raster; nodata handling; transforms.
# - geopandas: read AOI vector; CRS management; reprojection.
# - pandas: tabular cleaning, joins, Excel I/O helper.
# - numpy: fast array ops (unique SMUs, masking).
# - matplotlib (Agg backend): save QC plots in headless runs/CI.
# - shapely.mapping: convert AOI geometries to GeoJSON-like shapes for rasterio.mask.
# - openpyxl.Workbook: write multi-sheet Excel workbooks (Topsoil/Subsoil).



# -----------------------------
# User inputs (edit these only)
# -----------------------------
WORKING_DIR = Path(r"C:\Users\Asgharinia\OneDrive - Food and Agriculture Organization\Documents\PyAEZ\code\Soil_constraints")
RASTER_PATH  = WORKING_DIR / "HWSD2_RASTER/HWSD2.bil"
EXCEL_PATH   = WORKING_DIR / "HWSD2_LAYERS.xlsx"     # sheet: HWSD2_LAYERS
SHAPE_PATH   = WORKING_DIR / "Ghana/Ghana.shp"

OUT_TOPSOIL  = WORKING_DIR / "Top_Soil_Layers.xlsx"
OUT_SUBSOIL  = WORKING_DIR / "Sub_Soil_Layers.xlsx"
QC_RASTER_AOI_PNG = WORKING_DIR / "Raster_and_AOI_Boundary_Check.png"
QC_CLIP_PNG       = WORKING_DIR / "Clipped_Raster_Check.png"


# -----------------------------
# Verify inputs exist
# -----------------------------
for p in (RASTER_PATH, EXCEL_PATH, SHAPE_PATH):
    if not p.exists():
        raise FileNotFoundError(f"❌ File missing: {p}")
    else:
        print(f"✅ File found: {p}")


# -----------------------------
# Step 1: Load raster + AOI
# -----------------------------
with rasterio.open(RASTER_PATH) as src:
    raster_crs = src.crs
    raster_nodata = src.nodata if src.nodata is not None else -9999  # safer fallback
    print("Raster loaded successfully with CRS:", raster_crs)

# Load AOI; require CRS and reproject to raster CRS
aoi = gpd.read_file(SHAPE_PATH)
if aoi.crs is None:
    raise ValueError("AOI shapefile has no CRS defined. Please define its CRS before running.")
aoi = aoi.to_crs(raster_crs)

# Optional QC plot (uses a low-res read to avoid RAM hit)
with rasterio.open(RASTER_PATH) as src:
    # small, downsampled preview for plotting
    preview = src.read(1, out_shape=(1, max(1, src.height // 8), max(1, src.width // 8)))
plt.figure()
plt.imshow(preview, cmap="viridis")
aoi.boundary.plot(edgecolor="red", linewidth=2, ax=plt.gca())
plt.title("Raster and AOI Boundary Check")
plt.axis("off")
plt.savefig(QC_RASTER_AOI_PNG, dpi=300, bbox_inches="tight")
plt.close()


# -----------------------------------------
# Step 2: Mask raster with AOI polygon
# -----------------------------------------
with rasterio.open(RASTER_PATH) as src:
    shapes = [mapping(geom) for geom in aoi.geometry]
    clipped, clipped_transform = mask(src, shapes, crop=True, nodata=raster_nodata)
    clipped = clipped[0]
print("Raster clipped to AOI boundary.")

# Optional QC plot of clipped raster
plt.figure()
plt.imshow(clipped, cmap="viridis")
plt.title("Clipped Raster Check")
plt.colorbar()
plt.axis("off")
plt.savefig(QC_CLIP_PNG, dpi=300, bbox_inches="tight")
plt.close()


# -----------------------------------------
# Step 3: Extract unique SMUs inside AOI
# (fast path—no big dataframe needed)
# -----------------------------------------
valid = clipped != raster_nodata
unique_smu = np.unique(clipped[valid])
print("Unique soil mapping units extracted:", len(unique_smu))


# -----------------------------------------
# Step 4: Load & filter soil attributes
# -----------------------------------------
soil_attributes = pd.read_excel(EXCEL_PATH, sheet_name="HWSD2_LAYERS")
soil_data = soil_attributes[soil_attributes["HWSD2_SMU_ID"].isin(unique_smu)]
print("Soil attributes loaded and filtered by raster values.")


# -----------------------------------------
# Step 5: Clean & subset
# -----------------------------------------
soil_data = soil_data.fillna(0)
soil_data = soil_data[soil_data["TEXTURE_USDA"] != 0]
soil_data = soil_data.loc[soil_data.groupby(["HWSD2_SMU_ID", "LAYER"])["SHARE"].idxmax()]
print("Data cleaned and duplicates removed.")


# -----------------------------------------
# Step 6: Replace texture codes with labels
# -----------------------------------------
texture_lookup = pd.DataFrame({
    "CODE": range(1, 14),
    "VALUE": [
        "Clay (heavy)", "Silty clay", "Clay (light)", "Silty clay loam",
        "Clay loam", "Silt", "Silt loam", "Sandy clay", "Loam",
        "Sandy clay loam", "Sandy loam", "Loamy sand", "Sand"
    ]
})
soil_data = soil_data.merge(texture_lookup, left_on="TEXTURE_USDA", right_on="CODE", how="left")
soil_data["TEXTURE_USDA"] = soil_data["VALUE"]
soil_data.drop(columns=["CODE", "VALUE"], inplace=True)
print("Texture descriptions assigned.")


# -----------------------------------------
# Step 7: (OSD removed in this version)
# -----------------------------------------


# -----------------------------------------
# Step 8: Binary transform for ADD_PROP
# -----------------------------------------
soil_data["ADD_PROP"] = (soil_data["ADD_PROP"] == 3).astype(int)


# -----------------------------------------
# Step 9: Reorder & rename columns
# -----------------------------------------
column_order = [
    "HWSD2_SMU_ID", "TEXTURE_USDA", "ORG_CARBON", "PH_WATER", "TEB", "BSAT",
    "CEC_SOIL", "CEC_CLAY", "ROOT_DEPTH", "PHASE1", "PHASE2", "ROOTS", "IL",
    "DRAINAGE", "ESP", "ELEC_COND", "TCARBON_EQ", "GYPSUM", "COARSE",
    "ADD_PROP", "LAYER", "SHARE"
]
soil_data = soil_data[column_order]
soil_data.columns = [
    "CODE", "TXT", "OC", "pH", "TEB", "BS", "CEC_soil", "CEC_clay", "RSD",
    "SPR", "SPH", "ROOTS", "IL", "DRG", "ESP", "EC", "CCB", "GYP", "GRC",
    "VSP", "LAYER", "SHARE"
]
print("Columns reordered and renamed.")


# -----------------------------------------
# Step 10: Write Excel workbooks
# -----------------------------------------
def create_soil_workbook(data: pd.DataFrame, layers: list[str], filename: Path) -> None:
    wb = Workbook()
    wb.remove(wb.active)  # remove default empty sheet
    for layer in layers:
        layer_df = data[data["LAYER"] == layer]
        if not layer_df.empty:
            ws = wb.create_sheet(title=layer)
            ws.append(list(layer_df.columns))
            for row in layer_df.itertuples(index=False):
                ws.append(list(row))
    wb.save(filename)

top_layers = ["D1", "D2", "D3"]
sub_layers = ["D4", "D5", "D6", "D7"]
create_soil_workbook(soil_data, top_layers, OUT_TOPSOIL)
create_soil_workbook(soil_data, sub_layers, OUT_SUBSOIL)

print(f"✅ Wrote: {OUT_TOPSOIL.name} and {OUT_SUBSOIL.name}")

