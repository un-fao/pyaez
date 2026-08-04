# -*- coding: utf-8 -*-
"""
PyAEZ v3.0 (2026)
Module 4 — Step 1: Extract HWSD soil layers and slope distribution for selected AOI

Purpose
-------
This script:
    
1) HWSD: Downloads/prepares the original HWSD2 SMU raster from FAO link.
2) AOI: Downloads/prepares GAUL 2024 L2 boundary data from FAO link.
3) Selects an AOI using a user-defined boundary attribute and AOI value.
4) Clips the HWSD SMU raster to the AOI.
5) Extracts unique SMUs inside the AOI.
6) Filters HWSD2 layer attributes for those SMUs.
7) Reads the 10 GAEZ slope-share rasters:
      slp_cl01_30s ... slp_cl10_30s
8) Aligns the 10 slope-share rasters to the clipped HWSD AOI grid.
9) Reads sluslope30_hwsd2v10.dat to link:
      SMU × SEQ × slope class
10) Builds the final soil_slope_AOI table with:
      HWSD soil attributes
      slope class
      slope class definition
      AOI slope share
      SMU-SEQ-slope distribution weight
11) Exports one Excel workbook with:
      - Top_Soil_Layers
      - Sub_Soil_Layers
      - Slope_Class_Definition
      - SMU_Slope_Distribution
12) Exports key AOI rasters:
      - AOI HWSD SMU raster
      - AOI dominant slope-class raster

Important
---------
This version does NOT use GAEZ-V5.SLOPE-MED.tif.
The slope association is derived from the 10 slope-share rasters and
sluslope30_hwsd2v10.dat.

Optimized copy (v3_opt)
-----------------------
Same science as PyAEZ_Module4_Step1_Soil_Slope_v3.py; faster I/O and
vectorized aggregations only. See OPTIMIZATION_NOTES.md. Key opts marked
below with "# OPT:" comments — do not change weight/resample formulas.

Author / Credits
-----------------------------------------------------------------------------
Development and technical coordination:
 Shahla Asgharinia
 Food and Agriculture Organization of the United Nations (FAO)
 Land and Water Division, Geospatial Unit (NSLD)

Performance optimization:
 Mississippi State University
-----------------------------------------------------------------------------
"""
# -----------------------------
# -----------------------------
# %% Imports
# -----------------------------

from __future__ import annotations

from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import shutil
from pathlib import Path
import zipfile
import warnings
import time

import rasterio
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling

import geopandas as gpd
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors

from shapely.geometry import mapping, box


# -----------------------------
# -----------------------------
# %%  Runtime timer
# -----------------------------

SCRIPT_START_TIME = time.perf_counter()


def _stage(label: str, t0=None):
    # OPT: lightweight stage timer — call with t0=None to start, then again
    # with the returned clock to print elapsed seconds for that stage.
    """Print stage timing when PROFILE_STAGES is True. Returns new t0."""
    now = time.perf_counter()
    if PROFILE_STAGES and t0 is not None:
        print(f"[PROFILE] {label}: {now - t0:.2f}s")
    return now



# -----------------------------
# %%  USER INPUTS — EDIT THIS SECTION 
# Working directory
# -----------------------------
# OPT: FAO v3 hardcodes Shahla's OneDrive + "Inputs/". Resolve relative to
# this script so the Drive package ("input data/") works on any machine.
BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "input data"
OUTPUT_DIR = BASE_DIR / "Outputs" / "optimized"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# OPT: toggles for profiling / skipping matplotlib QC / slope NPZ cache.
# PROFILE_STAGES: print [PROFILE] wall times for major stages.
# SKIP_QC_PLOTS: skip PNG only (Excel/TIFs still written).
# USE_SLOPE_CACHE: reuse AOI-aligned 10-band slope stack on re-runs.
PROFILE_STAGES = True
SKIP_QC_PLOTS = True
USE_SLOPE_CACHE = True


# -----------------------------
# User inputs
# -----------------------------
# Select AOI using any valid field in the GAUL file.
# For country-level selection, use gaul0_name.
# For admin-1 selection, use gaul1_name.
# For admin-2 selection, use gaul2_name.
AOI_VALUES = ["United Republic of Tanzania"] # Area of interest, exact name used to select AOI from GAUL0
AOI_KEY = "tanzania"  # short name used for output files labeling



# %%  INPUTS

# OPT: Drive ships HWSD2_LAYERS_2.csv (~56 MB). Prefer CSV over forcing a
# slow xlsx round-trip; fall back to Excel if CSV is missing.
LAYERS_CSV_PATH = INPUT_DIR / "HWSD2_LAYERS_2.csv"
EXCEL_PATH = OUTPUT_DIR / "HWSD2_LAYERS_2.xlsx"
LAYERS_TABLE_PATH = LAYERS_CSV_PATH if LAYERS_CSV_PATH.exists() else EXCEL_PATH


#  Original HWSD2 SMU raster ZIP
HWSD_RASTER_URL = ("https://data.apps.fao.org/catalog/dataset/fedae0d8-677c-47a7-9049-6cba560e6c0a/resource/0e9b5791-2755-415f-8c26-b846a87e704e/download/hwsd2_raster.zip")
HWSD_ZIP = INPUT_DIR / "HWSD" / "hwsd2_raster.zip"
HWSD_DIR = INPUT_DIR / "HWSD" / "hwsd2_raster"


# AOI boundary source: GAUL 2024 L2
AOI_URL = "https://storage.googleapis.com/fao-maps-catalog-data/boundaries/GAUL_2024_L2.zip"

AOI_ZIP = INPUT_DIR / "Boundaries" / "GAUL_2024_L2.zip"
AOI_DIR = INPUT_DIR / "Boundaries" / "GAUL_2024_L2"
AOI_FIELD = "gaul0_name"
AOI_LABEL = AOI_KEY.lower().strip()


# # =============================================================================
# Slope-share inputs — GAEZ-style slope logic
# =============================================================================

SLOPE_SHARE_DIR = INPUT_DIR / "slp_cl"
SLUSLOPE_DAT = SLOPE_SHARE_DIR / "sluslope30_hwsd2v10.dat"

# =============================================================================
# HWSD / PyAEZ master grouping key
# =============================================================================
# SMU = AOI HWSD raster value
# COV, SEQ, SEQ2 = HWSD soil/component identity
# SLOPE_CLASS = GAEZ slope class
# LAYER = HWSD soil depth layer

Code = ["SMU", "COV", "SEQ", "SEQ2", "SLOPE_CLASS", "LAYER"]

TRACE_COLS = ["ID", "SMU", "SMU2", "SMU1", "COV", "SEQ", "SEQ2"]

DROP_COLS = ["Unnamed: 0"]   # keep ID
# -----------------------------
# %% Outputs
# -----------------------------
# Main tabular output
OUT_SOIL_SLOPE_AOI = OUTPUT_DIR / f"{AOI_LABEL}_soil_slope.xlsx"

# Key raster outputs
OUT_HWSD_SMU_TIF = OUTPUT_DIR / f"{AOI_LABEL}_hwsd_smu.tif"
OUT_SLOPE_CLASS_TIF = OUTPUT_DIR / f"{AOI_LABEL}_slope_class.tif"
OUT_SMU_SLOPE_CLASS_TIF = OUTPUT_DIR / f"{AOI_LABEL}_smu_slope_class.tif"  # multi band containing SMU-SEQ-SLP class
# QC PNG outputs
QC_RASTER_AOI_PNG = OUTPUT_DIR / f"raster_and_{AOI_LABEL}_boundary_check.png"
QC_CLIP_PNG = OUTPUT_DIR / f"clipped_hwsd_{AOI_LABEL}_check.png"
QC_SLOPE_PNG = OUTPUT_DIR / f"{AOI_LABEL}_slope_class_check.png"


# -----------------------------
# %% Helper functions
# -----------------------------


def download_if_needed(url: str, out_path: Path, label: str = "file") -> Path:
    """
    Download a file only if it does not already exist locally.
    Uses a browser-like User-Agent to avoid 403 errors from some servers.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"Using existing {label}: {out_path}")
        return out_path

    print(f"Downloading {label} from:\n{url}")

    try:
        request = Request(
            url,
            headers={
                "User-Agent": (
                    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0 Safari/537.36"
                )
            },
        )

        with urlopen(request, timeout=300) as response:
            with open(out_path, "wb") as f:
                shutil.copyfileobj(response, f)

        print(f"Downloaded {label} to: {out_path}")
        return out_path

    except HTTPError as e:
        raise RuntimeError(
            f"HTTP error while downloading {label}: {e.code} {e.reason}\n"
            f"URL:\n{url}\n"
            f"Target path:\n{out_path}"
        ) from e

    except URLError as e:
        raise RuntimeError(
            f"Network error while downloading {label}: {e.reason}\n"
            f"URL:\n{url}"
        ) from e
        
        
        


def unzip_if_needed(zip_path: Path, out_dir: Path, label: str = "archive") -> Path:
    """
    Extract a ZIP archive only if the output directory does not already contain files.
    Uses recursive search later, so it is robust to ZIPs that contain subfolders.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    existing_files = [p for p in out_dir.rglob("*") if p.is_file()]

    if not existing_files:
        print(f"Extracting {label}:\n{zip_path}")
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(out_dir)
        print(f"Extracted {label} to: {out_dir}")
    else:
        print(f"Using existing extracted {label}: {out_dir}")

    return out_dir


def find_raster_file(folder: Path, label: str = "raster") -> Path:
    """
    Find the main raster file inside an extracted raster folder.
    Priority:
      1) .bil files
      2) .tif / .tiff files

    This is suitable for the HWSD raster ZIP, which may contain BIL + sidecar files.
    """
    bil_files = sorted(folder.rglob("*.bil"))
    tif_files = sorted(folder.rglob("*.tif")) + sorted(folder.rglob("*.tiff"))

    candidates = bil_files or tif_files

    if not candidates:
        raise FileNotFoundError(
            f"No raster file (.bil, .tif, .tiff) found in extracted {label} folder:\n{folder}"
        )

    raster_path = candidates[0]
    print(f"Using {label}: {raster_path}")
    return raster_path


def find_vector_file(folder: Path, label: str = "boundary") -> Path:
    """
    Find the main vector file inside an extracted boundary folder.
    Priority:
      1) .shp files
      2) .gpkg files
      3) .geojson files
    """
    shp_files = sorted(folder.rglob("*.shp"))
    gpkg_files = sorted(folder.rglob("*.gpkg"))
    geojson_files = sorted(folder.rglob("*.geojson"))

    candidates = shp_files or gpkg_files or geojson_files

    if not candidates:
        raise FileNotFoundError(
            f"No vector file (.shp, .gpkg, .geojson) found in extracted {label} folder:\n{folder}"
        )

    vector_path = candidates[0]
    print(f"Using {label}: {vector_path}")
    return vector_path


def create_soil_slope_aoi_workbook(data: pd.DataFrame, filename: Path) -> None:
    """
    Create one Excel workbook containing two sheets:
      - Top_Soil_Layers: D1 only
      - Sub_Soil_Layers: D2 to D7
    """
    data = data.copy()
    data["LAYER"] = data["LAYER"].astype(str).str.upper().str.strip()

    topsoil_df = data[data["LAYER"] == "D1"].copy()
    subsoil_df = data[data["LAYER"].isin(["D2", "D3", "D4", "D5", "D6", "D7"])].copy()

    filename.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        topsoil_df.to_excel(
            writer,
            sheet_name="Top_Soil_Layers",
            index=False
        )

        subsoil_df.to_excel(
            writer,
            sheet_name="Sub_Soil_Layers",
            index=False
        )

    print(f"✅ Wrote single workbook: {filename.name}")
    print(f"   Sheet Top_Soil_Layers: {len(topsoil_df):,} rows")
    print(f"   Sheet Sub_Soil_Layers: {len(subsoil_df):,} rows")


def verify_paths(paths: list[Path]) -> None:
    """
    Verify required local paths exist.
    """
    for p in paths:
        if not Path(p).exists():
            raise FileNotFoundError(f"❌ Missing required input: {p}")
        print(f"✅ Found: {p}")

def standardize_hwsd_identity(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize HWSD identity fields.
    Step 1 keeps all AOI records.
    Step 2 will decide rating, zero suitability, or layer exclusion logic.
    """
    df = df.copy()

    # Harmonise SQ2 / SEQ2
    if "SEQ2" not in df.columns and "SQ2" in df.columns:
        df["SEQ2"] = df["SQ2"]

    # Convert identity columns where present
    id_cols = ["ID", "SMU2", "SMU1", "COV", "SEQ", "SEQ2", "SLOPE_CLASS"]

    for c in id_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

    # SMU is the AOI HWSD raster value; in this workflow it equals SMU2
    if "SMU2" in df.columns:
        df["SMU"] = df["SMU2"]

    return df


def add_master_code(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add one master Code column:
    SMU-COV-SEQ-SEQ2-SLOPE_CLASS-LAYER
    """
    df = df.copy()

    missing = [c for c in Code if c not in df.columns]

    if missing:
        raise ValueError(
            f"Cannot create CODE. Missing columns: {missing}"
        )

    df["CODE"] = (
        df[Code]
        .astype("string")
        .fillna("NA")
        .agg("-".join, axis=1)
    )

    return df
# -----------------------------
# %% Download / prepare inputs
# -----------------------------
download_if_needed(HWSD_RASTER_URL, HWSD_ZIP, label="original HWSD2 SMU raster ZIP")
unzip_if_needed(HWSD_ZIP, HWSD_DIR, label="original HWSD2 SMU raster")
RASTER_PATH = find_raster_file(HWSD_DIR, label="HWSD2 SMU raster")

download_if_needed(AOI_URL, AOI_ZIP, label="GAUL 2024 L2 boundary ZIP")
unzip_if_needed(AOI_ZIP, AOI_DIR, label="GAUL 2024 L2 boundaries")
SHAPE_PATH = find_vector_file(AOI_DIR, label="GAUL 2024 L2 boundary")



# -----------------------------
# %%  Verify required inputs
# -----------------------------
verify_paths([
    RASTER_PATH,
    LAYERS_TABLE_PATH,
    SHAPE_PATH,
    SLUSLOPE_DAT,
])

# =============================================================================
# %%  Slope-share helper functions
# =============================================================================

def find_slope_share_raster(folder: Path, slope_class: int) -> Path:
    """
    Find one GAEZ slope-share raster:
      slp_cl01_30s ... slp_cl10_30s

    The files may be extensionless, with .rdc sidecar files.
    """

    k = int(slope_class)

    candidates = [
        folder / f"slp_cl{k:02d}_30s",
        folder / f"slp_cl{k:02d}_30s.rst",
        folder / f"slp_cl{k:02d}_30s.tif",
        folder / f"slp_cl{k:02d}_30s.tiff",
        folder / f"slp_cl{k:02d}_30s.bil",
    ]

    for p in candidates:
        if p.exists() and p.is_file():
            return p

    hits = [
        p for p in folder.rglob(f"slp_cl{k:02d}_30s*")
        if p.is_file() and p.suffix.lower() != ".rdc"
    ]

    if hits:
        return sorted(hits)[0]

    raise FileNotFoundError(
        f"Missing slope-share raster for class {k:02d} in:\n{folder}"
    )


def list_slope_share_rasters(folder: Path) -> dict[int, Path]:
    """
    Return paths for all 10 slope-share rasters.
    """

    out = {}

    for k in range(1, 11):
        out[k] = find_slope_share_raster(folder, k)
        print(f"✅ Slope-share class {k:02d}: {out[k].name}")

    return out


def read_sluslope_dat(dat_path: Path, smu_filter=None) -> pd.DataFrame:
    """
    Read sluslope30_hwsd2v10.dat into long format.

    Output columns:
      SMU, SEQ, SLOPE_CLASS, SLUSLP_RAW, SLUSLP_WEIGHT

    SLUSLP_WEIGHT is normalized within each SMU + SLOPE_CLASS,
    so the sum across SEQ becomes 1.
    """
    # OPT: filter SMUs while reading. v3 parsed the whole global .dat into a
    # DataFrame then sliced; early continue + set lookup skips unused SMUs and
    # avoids a huge intermediate table. Weight math below is unchanged.

    smu_set = None
    if smu_filter is not None:
        smu_set = set(pd.Series(smu_filter).dropna().astype(int).tolist())

    rows = []

    with open(dat_path, "r", encoding="latin-1", errors="ignore") as f:
        for line_no, line in enumerate(f, start=1):

            if not line.strip():
                continue

            parts = line.split()

            try:
                # Common whitespace-separated format:
                # SMU SEQ AUX val1 ... val10
                if len(parts) >= 13:
                    smu = int(float(parts[0]))
                    # OPT: skip non-AOI SMUs before expanding 10 slope values
                    if smu_set is not None and smu not in smu_set:
                        continue
                    seq = int(float(parts[1]))
                    values = [float(x) for x in parts[3:13]]

                # Fixed-width fallback:
                # Fortran: (i5,i3,4x,10f8.6)
                else:
                    smu = int(line[0:5])
                    if smu_set is not None and smu not in smu_set:
                        continue
                    seq = int(line[5:8])
                    values = [
                        float(line[12 + i * 8:12 + (i + 1) * 8])
                        for i in range(10)
                    ]

            except Exception as e:
                raise ValueError(
                    f"Could not parse sluslope line {line_no}:\n{line}"
                ) from e

            # OPT: tuples are cheaper to append than per-row dicts
            for k, value in enumerate(values, start=1):
                rows.append((smu, seq, k, value))

    if not rows:
        raise ValueError(f"No rows read from {dat_path}")

    df = pd.DataFrame(rows, columns=["SMU", "SEQ", "SLOPE_CLASS", "SLUSLP_RAW"])

    df["SLUSLP_RAW"] = pd.to_numeric(df["SLUSLP_RAW"], errors="coerce").fillna(0.0)
    df.loc[df["SLUSLP_RAW"] < 0, "SLUSLP_RAW"] = 0.0

    denom = df.groupby(["SMU", "SLOPE_CLASS"])["SLUSLP_RAW"].transform("sum")

    df["SLUSLP_WEIGHT"] = np.where(
        denom > 0,
        df["SLUSLP_RAW"] / denom,
        0.0
    )

    return df.sort_values(
        ["SMU", "SEQ", "SLOPE_CLASS"],
        kind="mergesort"
    ).reset_index(drop=True)


def align_10_slope_share_rasters(
    slope_paths: dict[int, Path],
    reference_shape,
    reference_transform,
    reference_crs,
    valid_hwsd
) -> np.ndarray:
    """
    Align the 10 slope-share rasters to the clipped HWSD AOI grid.

    Returns:
      slope_stack with shape = (10, rows, cols)

    Values are normalized per pixel so that the 10 slope shares sum to 1.
    """
    # OPT: each slp_cl*.rst is ~global (~1.8 GB). v3 reprojected from the
    # full band; we crop to the AOI bbox (+ pad), then nearest-resample onto
    # the clipped HWSD grid — same Resampling.nearest as v3, less I/O.
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds, Window

    rows, cols = reference_shape
    raw_stack = np.zeros((10, rows, cols), dtype="float64")

    # Destination AOI extent in the HWSD/reference CRS
    left = reference_transform.c
    top = reference_transform.f
    right = left + cols * reference_transform.a
    bottom = top + rows * reference_transform.e

    for k in range(1, 11):

        aligned = np.full((rows, cols), np.nan, dtype="float64")

        with rasterio.open(slope_paths[k]) as src:
            # Map AOI bounds into the slope raster CRS, then read that window
            src_left, src_bottom, src_right, src_top = transform_bounds(
                reference_crs, src.crs, left, bottom, right, top, densify_pts=21
            )
            pad = abs(src.transform.a) * 2.0  # 2-pixel pad for edge safety
            win = from_bounds(
                src_left - pad,
                src_bottom - pad,
                src_right + pad,
                src_top + pad,
                transform=src.transform,
            )
            win = win.intersection(Window(0, 0, src.width, src.height))
            if win.width <= 0 or win.height <= 0:
                raise ValueError(f"Empty source window for slope class {k}")

            data = src.read(1, window=win, boundless=False, masked=False)
            win_transform = src.window_transform(win)

            reproject(
                source=data.astype("float64", copy=False),
                destination=aligned,
                src_transform=win_transform,
                src_crs=src.crs,
                src_nodata=src.nodata,
                dst_transform=reference_transform,
                dst_crs=reference_crs,
                dst_width=cols,
                dst_height=rows,
                dst_nodata=np.nan,
                resampling=Resampling.nearest,  # must stay nearest (v3)
            )

        aligned[~valid_hwsd] = np.nan
        aligned[aligned < 0] = np.nan

        raw_stack[k - 1] = aligned

    arr = np.where(np.isfinite(raw_stack), raw_stack, 0.0)
    total = arr.sum(axis=0)

    slope_stack = np.full_like(arr, np.nan)
    valid = valid_hwsd & (total > 0)

    # Same per-pixel normalize as v3 (vectorized across bands)
    slope_stack[:, valid] = arr[:, valid] / total[valid]

    return slope_stack


def build_smu_slope_distribution(
    clipped_smu: np.ndarray,
    valid_hwsd: np.ndarray,
    slope_stack: np.ndarray
) -> pd.DataFrame:
    """
    Build AOI-specific slope-class distribution for each SMU.

    Output:
      SMU, SLOPE_CLASS_MODE, SLOPE_CLASS_MODE_LABEL, SLP_share_1..SLP_share_10
    """
    # OPT: v3 looped each SMU with a full-grid boolean mask (O(n_SMU × pixels)).
    # bincount(weights=...) sums slope shares by SMU id in one pass per band —
    # same shares/mode, much less Python. Analogous to a weighted groupby.

    smu_flat = clipped_smu[valid_hwsd].astype(np.int64, copy=False)
    if smu_flat.size == 0:
        return pd.DataFrame(
            columns=["SMU", "SLOPE_CLASS_MODE", "SLOPE_CLASS_MODE_LABEL"]
            + [f"SLP_share_{k}" for k in range(1, 11)]
        )

    max_smu = int(smu_flat.max())
    totals = np.zeros((10, max_smu + 1), dtype=np.float64)
    for k in range(10):
        vals = slope_stack[k][valid_hwsd]
        vals = np.where(np.isfinite(vals), vals, 0.0)
        # totals[k, smu_id] = sum of band-k shares over pixels with that SMU
        totals[k] = np.bincount(smu_flat, weights=vals, minlength=max_smu + 1)

    smu_values = np.unique(smu_flat)
    records = []
    for smu in smu_values:
        smu = int(smu)
        t = totals[:, smu]
        total = float(t.sum())
        if total <= 0:
            continue
        shares = t / total
        mode_class = int(np.argmax(shares) + 1)
        rec = {
            "SMU": smu,
            "SLOPE_CLASS_MODE": mode_class,
            "SLOPE_CLASS_MODE_LABEL": SLP_CAPTIONS[mode_class],
        }
        for k in range(1, 11):
            rec[f"SLP_share_{k}"] = float(shares[k - 1])
        records.append(rec)

    return pd.DataFrame(records).sort_values(
        "SMU",
        kind="mergesort"
    ).reset_index(drop=True)
# =============================================================================
# %%  Step 1 — Load HWSD raster and select AOI
# =============================================================================
with rasterio.open(RASTER_PATH) as src:
    raster_crs = src.crs
    raster_nodata = src.nodata if src.nodata is not None else -9999
    raster_bounds = src.bounds

print("HWSD raster loaded successfully.")
print("Raster CRS:", raster_crs)
print("Raster NoData:", raster_nodata)
print("Raster bounds:", raster_bounds)


# Read GAUL boundary (AOI only when possible)
_t = _stage("start_gaul_load")
# OPT: largest baseline cost (~36s). v3 loads the full world L2 shapefile
# then filters in Python. where= pushes the predicate into the GIS driver
# (predicate pushdown) so we mainly load Tanzania polygons. Fallback if
# the driver ignores where= or returns empty.
_aoi_clauses = [f"{AOI_FIELD} = '{v}'" for v in AOI_VALUES]
_aoi_where = " OR ".join(_aoi_clauses)
try:
    aoi_global = gpd.read_file(SHAPE_PATH, where=_aoi_where)
    if aoi_global.empty:
        raise ValueError("where filter returned empty")
except Exception as _gaul_err:
    print(f"GAUL where-filter unavailable ({_gaul_err}); loading full file.")
    aoi_global = gpd.read_file(SHAPE_PATH)
_t = _stage("load GAUL boundaries", _t)

if aoi_global.empty:
    raise ValueError(f"AOI boundary file is empty: {SHAPE_PATH}")

if aoi_global.crs is None:
    raise ValueError(
        "AOI boundary file has no CRS defined. Define CRS before running, usually EPSG:4326."
    )

if AOI_FIELD not in aoi_global.columns:
    raise ValueError(
        f"AOI field '{AOI_FIELD}' not found in boundary file.\n"
        f"Available fields are:\n{list(aoi_global.columns)}"
    )


# Select AOI by user-defined field/value
aoi_targets = [v.strip().casefold() for v in AOI_VALUES]

aoi = aoi_global[
    aoi_global[AOI_FIELD]
    .astype(str)
    .str.strip()
    .str.casefold()
    .isin(aoi_targets)
].copy()

if aoi.empty:
    available = (
        aoi_global[AOI_FIELD]
        .dropna()
        .astype(str)
        .sort_values()
        .unique()
    )

    raise ValueError(
        f"No AOI found for AOI_VALUES={AOI_VALUES} using AOI_FIELD='{AOI_FIELD}'.\n"
        f"Check spelling or use another AOI_FIELD.\n"
        f"Example available values:\n{available[:40]}"
    )

print("\nSelected AOI:")
print(aoi[AOI_FIELD].drop_duplicates().to_string(index=False))
print(f"AOI output label: {AOI_LABEL}")


# Fix invalid geometries
aoi["geometry"] = aoi.geometry.buffer(0)

# Remove empty/null geometries
aoi = aoi[aoi.geometry.notna() & ~aoi.geometry.is_empty].copy()

if aoi.empty:
    raise ValueError("AOI geometries are empty after geometry cleaning.")


# Dissolve selected features into one AOI geometry
aoi = aoi.dissolve().reset_index(drop=True)

# Reproject AOI to HWSD raster CRS
aoi = aoi.to_crs(raster_crs)

# Check AOI intersects HWSD raster
raster_box = gpd.GeoDataFrame(
    geometry=[box(*raster_bounds)],
    crs=raster_crs
)

if not aoi.intersects(raster_box.geometry.iloc[0]).any():
    raise ValueError(
        "Selected AOI does not intersect the HWSD raster extent. "
        "Check CRS, AOI value, or input raster."
    )

print("AOI extracted, dissolved, and reprojected successfully.")
print("AOI bounds:", aoi.total_bounds)


# QC plot: raster preview + AOI boundary
with rasterio.open(RASTER_PATH) as src:
    preview = src.read(
        1,
        out_shape=(1, max(1, src.height // 8), max(1, src.width // 8))
    )

    extent = [
        src.bounds.left,
        src.bounds.right,
        src.bounds.bottom,
        src.bounds.top,
    ]




# =============================================================================
# %%  Step 2 — Clip HWSD raster with AOI
# =============================================================================
with rasterio.open(RASTER_PATH) as src:
    shapes = [mapping(geom) for geom in aoi.geometry]

    clipped, clipped_transform = mask(
        src,
        shapes,
        crop=True,
        nodata=raster_nodata
    )

    clipped = clipped[0]

print("HWSD raster clipped to AOI boundary.")


# =============================================================================
# %% Step 3 — Extract unique SMUs inside AOI
# =============================================================================
valid_hwsd = (clipped != raster_nodata) & np.isfinite(clipped) & (clipped > 0)
unique_smu = np.unique(clipped[valid_hwsd]).astype(int)

if unique_smu.size == 0:
    raise ValueError("No valid SMU values found inside AOI after clipping.")

print("\nUnique SMUs inside AOI:")
print(f"  Count: {len(unique_smu):,}")
print(f"  Min:   {unique_smu.min()}")
print(f"  Max:   {unique_smu.max()}")
print(f"  First values: {unique_smu[:20]}")

# -------------------------------------------------------------------------
# Export AOI-clipped HWSD SMU raster
# -------------------------------------------------------------------------
# This is the main SMU raster needed later for pixel-level map outputs.
# Each valid pixel contains the original HWSD SMU ID.

hwsd_to_write = np.where(
    valid_hwsd,
    clipped,
    raster_nodata
).astype("int32")

with rasterio.open(
    OUT_HWSD_SMU_TIF,
    "w",
    driver="GTiff",
    height=clipped.shape[0],
    width=clipped.shape[1],
    count=1,
    dtype="int32",
    crs=raster_crs,
    transform=clipped_transform,
    nodata=int(raster_nodata),
    compress="lzw",
) as dst:
    dst.write(hwsd_to_write, 1)
    dst.set_band_description(1, "HWSD_SMU")

print(f"✅ AOI HWSD SMU raster written: {OUT_HWSD_SMU_TIF.name}")

# =============================================================================
# %%  Step 4 — Load and filter HWSD layer attributes
# =============================================================================
_t = _stage("start_layers_load")
# OPT: read CSV when present (faster than openpyxl on ~56 MB); Excel fallback
if str(LAYERS_TABLE_PATH).lower().endswith(".csv"):
    soil_attributes = pd.read_csv(LAYERS_TABLE_PATH, low_memory=False)
else:
    soil_attributes = pd.read_excel(LAYERS_TABLE_PATH)
_t = _stage("load HWSD layers table", _t)

if "SMU2" not in soil_attributes.columns:
    raise ValueError(
        f"Column 'SMU2' not found in HWSD Excel file.\n"
        f"Available columns are:\n{list(soil_attributes.columns)}"
    )

soil_attributes["SMU2"] = pd.to_numeric(
    soil_attributes["SMU2"],
    errors="coerce"
).astype("Int64")

soil_data = soil_attributes[soil_attributes["SMU2"].isin(unique_smu)].copy()

# OPT: early filter — drop non-AOI SMUs before renames/joins (less RAM downstream)
_t = _stage("filter HWSD layers to AOI SMUs", _t)
print("\nHWSD attribute filtering:")
print(f"  Original rows: {len(soil_attributes):,}")
print(f"  Filtered rows: {len(soil_data):,}")

if soil_data.empty:
    raise ValueError(
        "No HWSD attribute rows matched the SMU values extracted from the raster. "
        "Check that the raster is the original SMU raster and that the Excel table uses matching SMU2 codes."
    )

soil_data["AOI_NAME"] = AOI_KEY
soil_data = standardize_hwsd_identity(soil_data)

soil_data["LAYER"] = soil_data["LAYER"].astype(str).str.upper().str.strip()
soil_data["LAYER_NUM"] = pd.to_numeric(
    soil_data["LAYER"].str.extract(r"(\d+)")[0],
    errors="coerce"
).astype("Int64")


# =============================================================================
# %%  Step 5 — Clean numeric columns
# =============================================================================
if "TX_US" not in soil_data.columns:
    warnings.warn("Column 'TX_US' not found.", RuntimeWarning)



_num_cols = [
    "SHARE", "ORG_C", "PH_W", "TOT_N", "CN_R", "TEB", "BSAT",
    "CEC_S", "CEC_C", "CEC_E", "AL_S", "ECE", "ESP", "TC_EQ", "GYPS",
    "CFR", "SAND", "SILT", "CLAY", "SEQ2", "BULK", "R_BULK",
    "RDEP", "AWCFL", "AWCFC", "CEC_C3"
]

for c in _num_cols:
    if c in soil_data.columns:
        soil_data[c] = pd.to_numeric(soil_data[c], errors="coerce")

print("HWSD attribute data cleaned. NaNs preserved.")


# =============================================================================
# %%  Step 6 — Texture labels
# =============================================================================
TX13_NAMES = {
    1: "Clay (heavy)",
    2: "Silty clay",
    3: "Clay",
    4: "Silty clay loam",
    5: "Clay loam",
    6: "Silt",
    7: "Silt loam",
    8: "Sandy clay",
    9: "Loam",
    10: "Sandy clay loam",
    11: "Sandy loam",
    12: "Loamy sand",
    13: "Sand",
}

def tx_label(x):
    if pd.isna(x):
        return None
    try:
        return TX13_NAMES.get(int(x))
    except Exception:
        return None

if "TX_US" in soil_data.columns:
    soil_data["TX_NAME"] = soil_data["TX_US"].map(tx_label)
    print("Texture labels added as TX_NAME.")


# =============================================================================
# %% Step 7 — Binary transform for ADD_P 
# =============================================================================
if "ADD_P" in soil_data.columns:
    soil_data["ADD_P"] = pd.to_numeric(soil_data["ADD_P"], errors="coerce")
    soil_data["VSP"] = np.where(soil_data["ADD_P"].eq(3), 1, 0)
    soil_data["GSP"] = np.where(soil_data["ADD_P"].eq(2), 1, 0)

# =============================================================================
# %%  Step 8 — Rename HWSD headers to planned PyAEZ names
# =============================================================================
rename_new_to_final = {"SQ2": "SEQ2",
    "TX_US": "TXT",
    "ORG_C": "OC",
    "PH_W": "pH",
    "BSAT": "BS",
    "CEC_S": "CEC_soil",
    "CEC_C": "CEC_clay",
    "ECE": "EC",
    "TC_EQ": "CCB",
    "GYPS": "GYP",
    "CFR": "GRC",
    "ROOTS": "ROO",
    "DRNG": "DRG",
    "RDEP": "ROOT_DEPTH",
    "PHS1": "PHASE1",
    "PHS2": "PHASE2",
}

safe_map = {
    old: new
    for old, new in rename_new_to_final.items()
    if old in soil_data.columns and (new not in soil_data.columns or new == old)
}

soil_out = soil_data.copy()
soil_out.rename(columns=safe_map, inplace=True)

print("Renamed HWSD columns:")
print(safe_map)


# =============================================================================
# =============================================================================
# %% Step 9 — Slope-share extraction and SMU–SEQ–slope link
# =============================================================================
# This replaces the old SLOPE-MED workflow.
#
# Final meaning:
#   SOIL_SLOPE_AOI rows = SMU × SEQ × SLOPE_CLASS × LAYER
# =============================================================================

SLP_CAPTIONS = {
    1: "0–0.5%",
    2: "0.5–2%",
    3: "2–5%",
    4: "5–8%",
    5: "8–12%",
    6: "12–16%",
    7: "16–24%",
    8: "24–30%",
    9: "30–45%",
    10: ">45%",
}


# -------------------------------------------------------------------------
# 9.1 Locate 10 slope-share rasters
# -------------------------------------------------------------------------

print("\nReading 10 slope-share rasters...")

slope_paths = list_slope_share_rasters(SLOPE_SHARE_DIR)

print("\nAligning 10 slope-share rasters to clipped HWSD grid...")

_t = _stage("start_slope_align")
# OPT: cache the AOI-aligned 10-band stack. Align is already fairly fast with
# windowed reads; NPZ helps repeated runs of the same AOI (invalidate if
# AOI grid / valid mask changes — shape + valid_hwsd checks below).
_cache_dir = OUTPUT_DIR / "cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
_slope_cache = _cache_dir / f"{AOI_LABEL}_slope_stack.npz"
_use_cache = False
if USE_SLOPE_CACHE and _slope_cache.exists():
    _cached = np.load(_slope_cache)
    if (
        _cached["slope_stack"].shape[1:] == clipped.shape
        and np.array_equal(_cached["valid_hwsd"], valid_hwsd)
    ):
        slope_stack = _cached["slope_stack"]
        _use_cache = True
        print(f"Using cached AOI slope stack: {_slope_cache.name}")
if not _use_cache:
    slope_stack = align_10_slope_share_rasters(
        slope_paths=slope_paths,
        reference_shape=clipped.shape,
        reference_transform=clipped_transform,
        reference_crs=raster_crs,
        valid_hwsd=valid_hwsd
    )
    if USE_SLOPE_CACHE:
        np.savez_compressed(
            _slope_cache, slope_stack=slope_stack, valid_hwsd=valid_hwsd
        )
        print(f"Wrote slope cache: {_slope_cache.name}")
_t = _stage("align_10_slope_share_rasters (+cache)", _t)

slope_sum = np.nansum(slope_stack, axis=0)
valid_slope = valid_hwsd & (slope_sum > 0)

if valid_slope.sum() == 0:
    raise ValueError("No valid slope-share pixels found inside AOI.")

print("\nSlope-share check")
print("-----------------")
print(f"Valid HWSD pixels:   {int(valid_hwsd.sum()):,}")
print(f"Valid slope pixels:  {int(valid_slope.sum()):,}")
print(f"Max abs(sum - 1):    {np.nanmax(np.abs(slope_sum[valid_slope] - 1)):.8f}")


# -------------------------------------------------------------------------
# 9.3 Export AOI dominant slope-class raster
# -------------------------------------------------------------------------
# This is a map/QC product.
# The table below still keeps all slope classes.

dominant_slope = np.full(clipped.shape, -9999, dtype="int16")

dominant_slope[valid_slope] = (
    np.nanargmax(slope_stack[:, valid_slope], axis=0) + 1
).astype("int16")

with rasterio.open(
    OUT_SLOPE_CLASS_TIF,
    "w",
    driver="GTiff",
    height=dominant_slope.shape[0],
    width=dominant_slope.shape[1],
    count=1,
    dtype="int16",
    crs=raster_crs,
    transform=clipped_transform,
    nodata=-9999,
    compress="lzw",
) as dst:
    dst.write(dominant_slope, 1)
    dst.set_band_description(1, "DOMINANT_SLOPE_CLASS")

print(f"✅ AOI slope-class raster written: {OUT_SLOPE_CLASS_TIF.name}")


# -------------------------------------------------------------------------
# 9.4 Slope QC map
# -------------------------------------------------------------------------
# OPT: matplotlib QC is optional; default skip saves a bit of time/IO.
# Does not affect Excel or GeoTIFF outputs used by Step 2.

if not SKIP_QC_PLOTS:
    plot_slope = dominant_slope.astype("float64")
    plot_slope[plot_slope == -9999] = np.nan

    cmap = plt.get_cmap("tab10", 10)
    bounds = np.arange(0.5, 11.5, 1.0)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    plt.figure(figsize=(8, 8))
    im = plt.imshow(plot_slope, cmap=cmap, norm=norm)
    cbar = plt.colorbar(im, ticks=np.arange(1, 11))
    cbar.ax.set_yticklabels([SLP_CAPTIONS[i] for i in range(1, 11)])
    cbar.set_label("Dominant slope class")
    plt.title(f"AOI dominant slope class — {AOI_KEY}")
    plt.axis("off")
    plt.savefig(QC_SLOPE_PNG, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"✅ Slope QC map saved: {QC_SLOPE_PNG.name}")
else:
    print("SKIP_QC_PLOTS=True — slope QC PNG not written.")


# -------------------------------------------------------------------------
# 9.5 Build AOI slope distribution by SMU
# -------------------------------------------------------------------------

_t = _stage("start_smu_slope_dist")
SMU_SLOPE_DISTRIBUTION = build_smu_slope_distribution(
    clipped_smu=clipped,
    valid_hwsd=valid_hwsd,
    slope_stack=slope_stack
)
_t = _stage("build_smu_slope_distribution", _t)

if SMU_SLOPE_DISTRIBUTION.empty:
    raise ValueError("SMU slope distribution table is empty.")

print("\nSMU slope distribution created.")
print(f"Rows: {len(SMU_SLOPE_DISTRIBUTION):,}")


# -------------------------------------------------------------------------
# 9.6 Read sluslope30_hwsd2v10.dat
# -------------------------------------------------------------------------

_t = _stage("start_sluslope_dat")
# OPT: unique_smu → early filter inside read_sluslope_dat (see function)
SLUSLOPE_LONG = read_sluslope_dat(
    SLUSLOPE_DAT,
    smu_filter=unique_smu
)
_t = _stage("read_sluslope_dat", _t)

if SLUSLOPE_LONG.empty:
    raise ValueError("No sluslope records found for AOI SMUs.")

print("\nsluslope file read.")
print(f"Rows: {len(SLUSLOPE_LONG):,}")


# =============================================================================
# %% Step 10 — Build final SOIL_SLOPE_AOI table
# =============================================================================
# Final table structure:
#   SMU × SEQ × SLOPE_CLASS × LAYER
#
# Each HWSD soil layer is repeated for each valid slope class associated with
# that SMU and soil component.
# =============================================================================

# Convert SMU slope distribution to long format
slope_share_cols = [f"SLP_share_{k}" for k in range(1, 11)]

SMU_SLOPE_LONG = SMU_SLOPE_DISTRIBUTION.melt(
    id_vars=["SMU"],
    value_vars=slope_share_cols,
    var_name="SLOPE_CLASS",
    value_name="AOI_SLOPE_SHARE"
)

SMU_SLOPE_LONG["SLOPE_CLASS"] = (
    SMU_SLOPE_LONG["SLOPE_CLASS"]
    .str.extract(r"(\d+)")
    .astype(int)
)

SMU_SLOPE_LONG["SLOPE_CLASS_LABEL"] = SMU_SLOPE_LONG["SLOPE_CLASS"].map(SLP_CAPTIONS)


# Combine AOI slope share with sluslope SMU–SEQ slope weights
SMU_SEQ_SLOPE = SLUSLOPE_LONG.merge(
    SMU_SLOPE_LONG,
    on=["SMU", "SLOPE_CLASS"],
    how="left"
)

SMU_SEQ_SLOPE["AOI_SLOPE_SHARE"] = SMU_SEQ_SLOPE["AOI_SLOPE_SHARE"].fillna(0.0)

SMU_SEQ_SLOPE["AOI_SEQ_SLOPE_WEIGHT"] = (
    SMU_SEQ_SLOPE["AOI_SLOPE_SHARE"] *
    SMU_SEQ_SLOPE["SLUSLP_WEIGHT"]
)

# Keep only meaningful combinations
SMU_SEQ_SLOPE = SMU_SEQ_SLOPE[
    (SMU_SEQ_SLOPE["AOI_SLOPE_SHARE"] > 0) &
    (SMU_SEQ_SLOPE["SLUSLP_WEIGHT"] > 0)
].copy()

if SMU_SEQ_SLOPE.empty:
    raise ValueError("No positive SMU–SEQ–slope combinations found.")


# AOI-specific total SEQ weight within each SMU
SEQ_WEIGHT = (
    SMU_SEQ_SLOPE.groupby(["SMU", "SEQ"], as_index=False)
    .agg(AOI_SEQ_WEIGHT=("AOI_SEQ_SLOPE_WEIGHT", "sum"))
)

denom = SEQ_WEIGHT.groupby("SMU")["AOI_SEQ_WEIGHT"].transform("sum")

SEQ_WEIGHT["AOI_SEQ_WEIGHT_NORM"] = np.where(
    denom > 0,
    SEQ_WEIGHT["AOI_SEQ_WEIGHT"] / denom,
    0.0
)

SMU_SEQ_SLOPE = SMU_SEQ_SLOPE.merge(
    SEQ_WEIGHT,
    on=["SMU", "SEQ"],
    how="left"
)


# Merge slope fields into HWSD soil layers
SOIL_SLOPE_AOI = soil_out.merge(
    SMU_SEQ_SLOPE[
        [
            "SMU",
            "SEQ",
            "SLOPE_CLASS",
            "SLOPE_CLASS_LABEL",
            "AOI_SLOPE_SHARE",
            "SLUSLP_RAW",
            "SLUSLP_WEIGHT",
            "AOI_SEQ_SLOPE_WEIGHT",
            "AOI_SEQ_WEIGHT",
            "AOI_SEQ_WEIGHT_NORM",
        ]
    ],
    on=["SMU", "SEQ"],
    how="left"
)

SOIL_SLOPE_AOI = standardize_hwsd_identity(SOIL_SLOPE_AOI)
SOIL_SLOPE_AOI = add_master_code(SOIL_SLOPE_AOI)

missing_slope = SOIL_SLOPE_AOI["SLOPE_CLASS"].isna().sum()

if missing_slope > 0:
    warnings.warn(
        f"{missing_slope:,} soil-layer rows have no slope-class match. "
        "Check sluslope records and AOI SMUs.",
        RuntimeWarning
    )


# -------------------------------------------------------------------------
# Clean and order final columns
# -------------------------------------------------------------------------

drop_cols = [c for c in DROP_COLS if c in SOIL_SLOPE_AOI.columns]

if drop_cols:
    SOIL_SLOPE_AOI = SOIL_SLOPE_AOI.drop(columns=drop_cols)


preferred_front = [
"AOI_NAME",
    "CODE",
    "ID",
    "SMU",
    "SMU2",
    "SMU1",
    "COV",
    "SEQ",
    "SEQ2",
    "SLOPE_CLASS",
    "SLOPE_CLASS_LABEL",
    "AOI_SLOPE_SHARE",
    "SLUSLP_RAW",
    "SLUSLP_WEIGHT",
    "AOI_SEQ_SLOPE_WEIGHT",
    "AOI_SEQ_WEIGHT",
    "AOI_SEQ_WEIGHT_NORM",
    "LAYER",
    "LAYER_NUM",
    "SHARE",
    "TXT",
    "TX_NAME",
    "SAND",
    "SILT",
    "CLAY",
    "FAO90",
    "F90ID",
    "SWR",
    "OC",
    "pH",
    "TEB",
    "BS",
    "CEC_soil",
    "CEC_clay",
    "EC",
    "ESP",
    "CCB",
    "GYP",
    "GRC",
    "PHASE1",
    "PHASE2",
    "ROOT_DEPTH",
    "ROO",
    "IL",
    "DRG",
    "ADD_P",
    "VSP",
    "GSP",
    "AWCFL",
    "AWCFC",
]

front = [c for c in preferred_front if c in SOIL_SLOPE_AOI.columns]
rest = [c for c in SOIL_SLOPE_AOI.columns if c not in front]

SOIL_SLOPE_AOI = SOIL_SLOPE_AOI[front + rest].sort_values(
    ["SMU", "SEQ", "SLOPE_CLASS", "LAYER_NUM"],
    kind="mergesort"
).reset_index(drop=True)

print("\nFinal SOIL_SLOPE_AOI table created.")
print(f"Rows:        {len(SOIL_SLOPE_AOI):,}")
print(f"Unique SMU:  {SOIL_SLOPE_AOI['SMU'].nunique(dropna=True):,}")
print(f"Unique CODE: {SOIL_SLOPE_AOI['CODE'].nunique(dropna=True):,}")


# =============================================================================
# %% Step 11 — Export Excel workbook
# =============================================================================

SLOPE_CLASS_DEFINITION = pd.DataFrame({
    "SLOPE_CLASS": list(SLP_CAPTIONS.keys()),
    "SLOPE_CLASS_LABEL": list(SLP_CAPTIONS.values())
})

_t = _stage("start_excel_export")
with pd.ExcelWriter(OUT_SOIL_SLOPE_AOI, engine="openpyxl") as writer:

    SOIL_SLOPE_AOI[
        SOIL_SLOPE_AOI["LAYER"].astype(str).str.upper().eq("D1")
    ].to_excel(
        writer,
        sheet_name="Top_Soil_Layers",
        index=False
    )

    SOIL_SLOPE_AOI[
        SOIL_SLOPE_AOI["LAYER"].astype(str).str.upper().isin(
            ["D2", "D3", "D4", "D5", "D6", "D7"]
        )
    ].to_excel(
        writer,
        sheet_name="Sub_Soil_Layers",
        index=False
    )

    SLOPE_CLASS_DEFINITION.to_excel(
        writer,
        sheet_name="Slope_Class_Definition",
        index=False
    )

    SMU_SLOPE_DISTRIBUTION.to_excel(
        writer,
        sheet_name="SMU_Slope_Distribution",
        index=False
    )

print(f"✅ Final Excel output: {OUT_SOIL_SLOPE_AOI}")
_t = _stage("excel export", _t)


# =============================================================================
# %% Step 12 — Quick checks
# =============================================================================

print("\nQuick checks")
print("------------")
print("Sample rows:")
print(SOIL_SLOPE_AOI.head(3).to_string(index=False))

if {"SAND", "SILT", "CLAY"}.issubset(SOIL_SLOPE_AOI.columns):
    sums = (
        SOIL_SLOPE_AOI["SAND"].fillna(0) +
        SOIL_SLOPE_AOI["SILT"].fillna(0) +
        SOIL_SLOPE_AOI["CLAY"].fillna(0)
    )
    print("SAND+SILT+CLAY off-by >3 rows:", int((sums - 100).abs().gt(3).sum()))

print("Unique (SMU, SEQ):", SOIL_SLOPE_AOI.drop_duplicates(["SMU", "SEQ"]).shape[0])
print(
    "Unique (SMU, SEQ, SLOPE_CLASS):",
    SOIL_SLOPE_AOI.drop_duplicates(["SMU", "SEQ", "SLOPE_CLASS"]).shape[0]
)

tmp_slope = (
    SOIL_SLOPE_AOI[["SMU", "SLOPE_CLASS", "AOI_SLOPE_SHARE"]]
    .drop_duplicates()
    .groupby("SMU")["AOI_SLOPE_SHARE"]
    .sum()
)

print("SMU slope-share sum min/max:", round(float(tmp_slope.min()), 6), round(float(tmp_slope.max()), 6))

tmp_seq = (
    SOIL_SLOPE_AOI[["SMU", "SEQ", "AOI_SEQ_WEIGHT_NORM"]]
    .drop_duplicates()
    .groupby("SMU")["AOI_SEQ_WEIGHT_NORM"]
    .sum()
)

print("SMU SEQ-weight sum min/max:", round(float(tmp_seq.min()), 6), round(float(tmp_seq.max()), 6))

print("Topsoil rows:", int((SOIL_SLOPE_AOI["LAYER"] == "D1").sum()))
print(
    "Subsoil rows:",
    int(SOIL_SLOPE_AOI["LAYER"].isin(["D2", "D3", "D4", "D5", "D6", "D7"]).sum())
)

print("\nMain outputs")
print("------------")
print(f"Tabular output: {OUT_SOIL_SLOPE_AOI}")
print(f"HWSD SMU map:   {OUT_HWSD_SMU_TIF}")
print(f"Slope map:      {OUT_SLOPE_CLASS_TIF}")


# =============================================================================
# %% Total computational time
# =============================================================================

SCRIPT_END_TIME = time.perf_counter()
elapsed_seconds = SCRIPT_END_TIME - SCRIPT_START_TIME

print("\n" + "=" * 60)
print(f"Total computational time: {elapsed_seconds:.2f} seconds")
print(f"Total computational time: {elapsed_seconds / 60:.2f} minutes")
print("=" * 60)

# -----------------------------
