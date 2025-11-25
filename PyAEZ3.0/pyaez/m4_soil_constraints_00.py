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
   source: https://data.apps.fao.org/catalog/iso/ff5c613c-75bb-46a9-a162-bc728059b465
2) Implements AOI-based SMU extraction for seven soil layers (D1–D7).
3) extracts the slope class using global slop database for the AoI
   source: https://data.apps.fao.org/catalog//iso/44950f32-a84d-4784-ac50-1148beeb8597
    
    Inputs
    - HWSD2 raster (BIL + sidecars)
    - AOI boundary (shp file)
    - HWSD attribute workbook (Excel: HWSD2_LAYERS_2)
    - SOILFER.SLOPE-MED.tif (global slope database)
    Outputs
    - Cleaned per-layer tables for D1–D7 (e.g., Topsoil D1, Subsoil D2–D7)
    - Optional QC plots for global soil map
    - Optional QC plots for soil map for AOI 
    - Optional QC plots for slop classes of AOI 
    

"""

# -----------------------------
# Imports
# -----------------------------
from pathlib import Path
import os, sys
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
from rasterio.warp import reproject, Resampling
import warnings
from matplotlib import colors


# -----------------------------
# Working directory
# -----------------------------
try:
    WORKING_DIR = Path(__file__).resolve().parent
except NameError:
    WORKING_DIR = Path.cwd()

# If running in notebooks/Spyder, force your project path:
WORKING_DIR = Path(r"C:\Users\Asgharinia\OneDrive - Food and Agriculture Organization\Documents\PyAEZ\code\Soil_constraints")

# -----------------------------
# User inputs
# -----------------------------
RASTER_PATH = WORKING_DIR / "hwsd2_raster" / "HWSD2.bil"
EXCEL_PATH  = WORKING_DIR / "HWSD2_LAYERS_2.xlsx"   # new HWSD layer file
SHAPE_PATH  = WORKING_DIR / "Tanzania" / "Tanzania.shp"
SLOPE_TIF   = WORKING_DIR / "SOILFER.SLOPE-MED.tif"  # categorical slope classes 0..10

OUT_TOPSOIL       = WORKING_DIR / "Top_Soil_Layers.xlsx"
OUT_SUBSOIL       = WORKING_DIR / "Sub_Soil_Layers.xlsx"
QC_RASTER_AOI_PNG = WORKING_DIR / "Raster_and_AOI_Boundary_Check.png"
QC_CLIP_PNG       = WORKING_DIR / "Clipped_Raster_Check.png"
QC_SLOPE_PNG      = WORKING_DIR / "Clipped_Slope_Check.png"
OUT_SLP_CSV       = WORKING_DIR / "SMU_SLOPE_CLASSES.csv"

# -----------------------------
# Verify inputs
# -----------------------------
for p in (RASTER_PATH, EXCEL_PATH, SHAPE_PATH, SLOPE_TIF):
    if not p.exists():
        raise FileNotFoundError(f"❌ File missing: {p}")
    else:
        print(f"✅ File found: {p}")





# -----------------------------
# Step 1 – Load HWSD raster + AOI
# -----------------------------
with rasterio.open(RASTER_PATH) as src:
    raster_crs    = src.crs
    raster_nodata = src.nodata if src.nodata is not None else -9999
    print("Raster loaded successfully with CRS:", raster_crs)

# Load AOI; reproject to raster CRS
aoi = gpd.read_file(SHAPE_PATH)
if aoi.crs is None:
    raise ValueError("AOI shapefile has no CRS defined. Please define its CRS before running.")
aoi = aoi.to_crs(raster_crs)

# QC plot: raster preview + AOI boundary
with rasterio.open(RASTER_PATH) as src:
    preview = src.read(
        1,
        out_shape=(1, max(1, src.height // 8), max(1, src.width // 8))
    )

plt.figure()
plt.imshow(preview, cmap="viridis")
aoi.boundary.plot(edgecolor="red", linewidth=2, ax=plt.gca())
plt.title("Raster and AOI Boundary Check")
plt.axis("off")
plt.savefig(QC_RASTER_AOI_PNG, dpi=300, bbox_inches="tight")
plt.close()

# -----------------------------
# Step 2 – Mask HWSD raster with AOI
# -----------------------------
with rasterio.open(RASTER_PATH) as src:
    shapes = [mapping(geom) for geom in aoi.geometry]
    clipped, clipped_transform = mask(src, shapes, crop=True, nodata=raster_nodata)
    clipped = clipped[0]

print("Raster clipped to AOI boundary.")

plt.figure()
plt.imshow(clipped, cmap="viridis")
plt.title("Clipped Raster Check")
plt.colorbar()
plt.axis("off")
plt.savefig(QC_CLIP_PNG, dpi=300, bbox_inches="tight")
plt.close()

# -----------------------------
# Step 3 – Unique SMUs inside AOI
# -----------------------------
valid = clipped != raster_nodata
unique_smu = np.unique(clipped[valid])
print("Unique soil mapping units in AOI:", len(unique_smu))

# -----------------------------
# Step 4 – Load & filter HWSD layers 
# -----------------------------
soil_attributes = pd.read_excel(EXCEL_PATH)  # expects new headers listed above

# Filter by SMU2 (matches raster SMU ids)
soil_data = soil_attributes[soil_attributes["SMU2"].isin(unique_smu)].copy()
print("Soil attributes rows after SMU filter:", len(soil_data))

# Augment IDs
soil_data["SMU"] = pd.to_numeric(soil_data["SMU2"], errors="coerce").astype("Int64")
soil_data["SEQ"] = pd.to_numeric(soil_data["SEQ"], errors="coerce").astype("Int64")
soil_data["CODE"] = soil_data["SMU"].astype(str) + "-" + soil_data["SEQ"].astype(str)

soil_data["LAYER"] = soil_data["LAYER"].astype(str).str.upper().str.strip()
soil_data["LAYER_NUM"] = pd.to_numeric(
    soil_data["LAYER"].str.extract(r"(\d+)")[0],
    errors="coerce"
).astype("Int64")

# -----------------------------
# Step 5 – Clean & subset (minimal)
# -----------------------------
# keep rows with texture code
soil_data = soil_data[soil_data["TX_US"].notna()]

# Light numeric coercion; no NaN filling
_num_cols = [
    "SHARE","ORG_C","PH_W","TOT_N","CN_R","TEB","BSAT",
    "CEC_S","CEC_C","CEC_E","AL_S","ECE","ESP","TC_EQ","GYPS",
    "CFR","SAND","SILT","CLAY","SQ2","BULK","R_BULK",
    "RDEP","AWCFL","AWCFC","CEC_C3"
]
for c in _num_cols:
    if c in soil_data.columns:
        soil_data[c] = pd.to_numeric(soil_data[c], errors="coerce")

print("Data cleaned (no component collapse; NaNs preserved).")

# -----------------------------
# Step 6 – Texture labels (TX_US → TX_NAME)
# -----------------------------
TX13_NAMES = {
    1:"Clay (heavy)",  2:"Silty clay",    3:"Clay",
    4:"Silty clay loam", 5:"Clay loam",   6:"Silt",
    7:"Silt loam",     8:"Sandy clay",    9:"Loam",
    10:"Sandy clay loam", 11:"Sandy loam", 12:"Loamy sand", 13:"Sand"
}

def tx_label(x):
    if pd.isna(x):
        return None
    try:
        return TX13_NAMES.get(int(x))
    except Exception:
        return None

soil_data["TX_NAME"] = soil_data["TX_US"].map(tx_label)
print("Texture labels added as TX_NAME; TX_US kept as numeric code.")

# -----------------------------
# Step 7 – Binary transform for ADD_P (was ADD_PROP)
# -----------------------------
# 3 = presence; other = 0
if "ADD_P" in soil_data.columns:
    soil_data["ADD_P"] = (soil_data["ADD_P"] == 3).astype(int)

# -----------------------------
# Step 8 – Map new HWSD headers → planned names, keep ALL columns
# -----------------------------
# Map from new HWSD headers -> final planned headers
rename_new_to_final = {
    # texture
    "TX_US": "TXT",

    # chem/phys (short planned labels)
    "ORG_C": "OC",
    "PH_W": "pH",
    "BSAT": "BS",
    "CEC_S": "CEC_soil",
    "CEC_C": "CEC_clay",
    "ECE": "EC",
    "TC_EQ": "CCB",
    "GYPS": "GYP",
    "CFR": "GRC",

    # phases & classes
    "ROOTS": "ROO",
    "DRNG": "DRG",
    "ADD_P": "VSP",
    "RDEP": "ROOT_DEPTH",
    "PHS1": "PHASE1",
    "PHS2": "PHASE2",
}

# Safe renaming (don't overwrite existing targets)
safe_map = {
    old: new for old, new in rename_new_to_final.items()
    if old in soil_data.columns and (new not in soil_data.columns or new == old)
}

soil_out = soil_data.copy()
soil_out.rename(columns=safe_map, inplace=True)

# Preferred columns at front; keep everything else after
preferred_front = [
    # identity & weights
    "CODE","SMU","SEQ","LAYER","LAYER_NUM","SHARE",
    # texture
    "TXT","TX_NAME","SAND","SILT","CLAY","FAO90","SWR",
    # chem/phys
    "OC","pH","TEB","BS","CEC_soil","CEC_clay","EC","ESP","CCB","GYP","GRC",
    # phases & classes
    "PHASE1","PHASE2","ROOT_DEPTH","ROO","IL","DRG","VSP",
    # slope (will be added later)
    "SLP","SLP_caption","slope_pct_mean_est","slope_pct_median_est"
]

front = [c for c in preferred_front if c in soil_out.columns]
rest  = [c for c in soil_out.columns if c not in front]
soil_out = soil_out.loc[:, front + rest]

print("Renamed (new -> final):", safe_map)
missing_pref = [c for c in preferred_front if c not in soil_out.columns]
if missing_pref:
    warnings.warn(f"Preferred columns missing (kept all others): {missing_pref}", RuntimeWarning)

# -----------------------------
# Step 9 – Slope classes from SOILFER.SLOPE-MED (categorical)
# -----------------------------


# FAO/GAEZ 10-class scheme (class 1..10)
SLP_CAPTIONS = {
    1:"0–0.5%",  2:"0.5–2%",  3:"2–5%",  4:"5–8%",  5:"8–12%",
    6:"12–16%",  7:"16–24%",  8:"24–30%", 9:"30–45%", 10:">45%"
}
SLP_MIDPOINTS = np.array(
    [0.25, 1.25, 3.5, 6.5, 10.0, 14.0, 20.0, 27.0, 37.5, 50.0],
    dtype=float
)

# ---- 9.1 Clip slope raster by AOI polygon (same way as HWSD) ----
with rasterio.open(SLOPE_TIF) as slp_ds:
    slope_crs   = slp_ds.crs
    slope_nd    = slp_ds.nodata

    # AOI in slope CRS (in case it differs from HWSD CRS)
    aoi_slope = aoi.to_crs(slope_crs)

    slope_clip, slope_transform = mask(
        slp_ds,
        [mapping(g) for g in aoi_slope.geometry],
        crop=True,
        nodata=slope_nd
    )
    slope_clip = slope_clip[0].astype("float64")

# Treat explicit nodata as NaN
if slope_nd is not None:
    slope_clip[slope_clip == slope_nd] = np.nan

# GAEZ slope products use 1..10 for classes and 0 (or 255) as "no data".
# We therefore treat <=0 as nodata:
slope_clip[slope_clip <= 0] = np.nan

# ---- 9.2 Align slope clip to HWSD clip grid  ----
if slope_clip.shape == clipped.shape:
    # grids already match 
    slope_on_hwsd = slope_clip
else:
    slope_on_hwsd = np.full(clipped.shape, np.nan, dtype="float64")
    with rasterio.open(SLOPE_TIF) as slp_ds:
        reproject(
            source=slope_clip,
            destination=slope_on_hwsd,
            src_transform=slope_transform,
            src_crs=slp_ds.crs,
            dst_transform=clipped_transform,
            dst_crs=raster_crs,
            dst_width=clipped.shape[1],
            dst_height=clipped.shape[0],
            resampling=Resampling.nearest,
            src_nodata=np.nan,
            dst_nodata=np.nan
        )

# Restrict strictly to HWSD valid pixels inside AOI
slp_class = slope_on_hwsd.copy()
slp_class[clipped == raster_nodata] = np.nan   # outside AOI/HWSD

# At this point slp_class contains integer codes 1..10 (NaN outside AOI or nodata)

# QC plot: show real slope ranges in legend
cmap  = plt.get_cmap("tab10", 10)                       # 10 discrete colors
bounds = np.arange(0.5, 11.5, 1.0)                      # np.arange(start, stop, step), class boundaries: 0.5,1.5,...,10.5
norm   = colors.BoundaryNorm(bounds, cmap.N)

plt.figure(figsize=(8, 8))
im = plt.imshow(slp_class, cmap=cmap, norm=norm)
plt.title("Slope (%) classes (SOILFER.SLOPE-MED) – AOI clipped")
cbar = plt.colorbar(im, ticks=np.arange(1, 11))

# Use the real ranges as labels
cbar.ax.set_yticklabels([SLP_CAPTIONS[i] for i in range(1, 11)])
cbar.set_label("Slope (%)")

plt.axis("off")
plt.savefig(QC_SLOPE_PNG, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved slope QC to: {QC_SLOPE_PNG.name}")


# ---- 9.3 Aggregate per SMU (mode, shares, mean/median %) ----
valid_slope = np.isfinite(slp_class) & (clipped != raster_nodata)

smu_vals = clipped[valid_slope].astype(int)
cls_vals = slp_class[valid_slope].astype(int)

records = []
for smu in np.unique(smu_vals):
    mask_smu = smu_vals == smu
    classes_smu = cls_vals[mask_smu]
    if classes_smu.size == 0:
        continue

    # counts for classes 1..10 (index 0 unused)
    counts = np.bincount(classes_smu, minlength=11)
    if counts[1:].sum() == 0:
        continue

    # dominant slope class for this SMU
    mode_class = int(np.argmax(counts[1:]) + 1)

    rec = {
        "SMU": int(smu),
        "SLP": mode_class,
        "SLP_caption": SLP_CAPTIONS.get(mode_class)
    }

    total = counts[1:].sum()

    # fractional area share per class (optional but nice to keep)
    for k in range(1, 11):
        rec[f"SLP_share_{k}"] = float(counts[k] / total)

    # approximate mean/median slope % for this SMU using class midpoints
    mid_vals = SLP_MIDPOINTS[classes_smu - 1]
    rec["slope_pct_mean_est"]   = float(np.nanmean(mid_vals))
    rec["slope_pct_median_est"] = float(np.nanmedian(mid_vals))

    records.append(rec)

SLP_BY_SMU = (
    pd.DataFrame.from_records(records)
      .sort_values("SMU")
      .reset_index(drop=True)
)

SLP_BY_SMU.to_csv(OUT_SLP_CSV, index=False)
print(f"✅ Wrote slope-by-SMU table: {OUT_SLP_CSV.name} (rows: {len(SLP_BY_SMU)})")

# ---- 9.4 Merge SLP into soil_out ----
soil_out = soil_out.merge(
    SLP_BY_SMU[["SMU", "SLP", "SLP_caption",
                "slope_pct_mean_est", "slope_pct_median_est"]],
    on="SMU", how="left"
)



# -----------------------------
# Final cleanup: slope + IDs + SMU columns
# -----------------------------


# 1) Drop stray index / ID columns
drop_cols = [c for c in ["Unnamed: 0", "ID"] if c in soil_out.columns]
if drop_cols:
    soil_out = soil_out.drop(columns=drop_cols)

# 2) Clean up duplicated slope columns from merge
#    (keep x, drop y, and rename x to final names)
cols_to_drop = [c for c in ["SLP_y", "SLP_caption_y",
                            "slope_pct_mean_est_y", "slope_pct_median_est_y"]
                if c in soil_out.columns]
if cols_to_drop:
    soil_out = soil_out.drop(columns=cols_to_drop)

rename_map = {}
if "SLP_x" in soil_out.columns:
    rename_map["SLP_x"] = "SLOPE_CLASS"       # your requested final name
if "SLP_caption_x" in soil_out.columns:
    rename_map["SLP_caption_x"] = "SLP_caption"
if "slope_pct_mean_est_x" in soil_out.columns:
    rename_map["slope_pct_mean_est_x"] = "slope_pct_mean_est"
if "slope_pct_median_est_x" in soil_out.columns:
    rename_map["slope_pct_median_est_x"] = "slope_pct_median_est"

if rename_map:
    soil_out = soil_out.rename(columns=rename_map)

# 3) Consolidate SMU / SMU1 / SMU2
if {"SMU", "SMU1", "SMU2"}.issubset(soil_out.columns):
    # make sure they are comparable
    for col in ["SMU", "SMU1", "SMU2"]:
        soil_out[col] = pd.to_numeric(soil_out[col], errors="coerce").astype("Int64")

    same_all = (soil_out["SMU"].eq(soil_out["SMU1"]) &
                soil_out["SMU"].eq(soil_out["SMU2"]))

    if same_all.all():
        # all three columns are identical everywhere -> keep single SMU
        soil_out = soil_out.drop(columns=["SMU1", "SMU2"])
        print("SMU, SMU1, SMU2 identical for all rows → keeping only SMU.")
    else:
        # some differences exist → keep as-is so you can inspect
        n_diff = (~same_all).sum()
        print(f"SMU, SMU1, SMU2 differ in {n_diff} rows → keeping all 3 columns.")


# -----------------------------
# Step 10 – Write Excel workbooks (Topsoil/Subsoil)
# -----------------------------
def create_soil_workbook(data: pd.DataFrame, layers: list[str], filename: Path) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    for layer in layers:
        layer_df = data[data["LAYER"] == layer]
        if not layer_df.empty:
            ws = wb.create_sheet(title=layer)
            ws.append(list(layer_df.columns))
            for row in layer_df.itertuples(index=False):
                ws.append(list(row))
    wb.save(filename)

top_layers = ["D1"]
sub_layers = ["D2","D3","D4","D5","D6","D7"]

create_soil_workbook(soil_out, top_layers, OUT_TOPSOIL)
create_soil_workbook(soil_out, sub_layers, OUT_SUBSOIL)

print(f"✅ Wrote: {OUT_TOPSOIL.name} and {OUT_SUBSOIL.name}")

# -----------------------------
# Quick checks
# -----------------------------
print("Sample rows:\n", soil_out.head(3))
if {"SAND","SILT","CLAY"}.issubset(soil_out.columns):
    sums = (soil_out["SAND"].fillna(0) +
            soil_out["SILT"].fillna(0) +
            soil_out["CLAY"].fillna(0))
    print("SAND+SILT+CLAY off-by>3 in rows:", (sums-100).abs().gt(3).sum())

print("Unique (SMU,SEQ):", soil_out.drop_duplicates(["SMU","SEQ"]).shape[0])
print("Done.")
