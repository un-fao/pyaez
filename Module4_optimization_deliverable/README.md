# PyAEZ Module 4 Step 1 — Optimization Package

Faster AOI soil/slope table builds with **the same science** as FAO `PyAEZ_Module4_Step1_Soil_Slope_v3.py`.  
Validated: Tanzania — Excel keys/attributes match; weight columns within ~1e-14 float noise; SMU and slope GeoTIFFs bit-identical.  
Runtime (reference PC): ~65 s → ~22 s (~2.9×).

**Scope:** Step 1 only. Step 2 is unchanged (still FAO’s unfinished v18).

## Package contents

| Path | Purpose |
|------|---------|
| `PyAEZ_Module4_Step1_Soil_Slope_v3_opt.py` | Optimized Step 1 script (drop-in alongside original v3) |
| `docs/OPTIMIZATION_NOTES.md` | Technical summary (timings, changes, principles) |
| `docs/OPTIMIZATION_NOTES_SIMPLE.md` | Non-technical summary for reports / email |

## Install

1. Keep the original FAO script for comparison:
   - `PyAEZ_Module4_Step1_Soil_Slope_v3.py` (do not overwrite)
2. Copy `PyAEZ_Module4_Step1_Soil_Slope_v3_opt.py` into your `PyAEZ3.0/module4/` folder.
3. Point inputs at your Drive package (or keep the script’s defaults):
   - Script folder as `BASE_DIR`
   - `input data/` for HWSD, GAUL, `slp_cl/`, `HWSD2_LAYERS_2.csv`
4. Run with the `pyaez` conda env (needs `rasterio`, `geopandas`, `pandas`, `openpyxl`).

```text
conda activate pyaez
python PyAEZ_Module4_Step1_Soil_Slope_v3_opt.py
```

Outputs go to `Outputs/optimized/` by default.

## Optional flags (top of `_opt.py`)

| Flag | Default | Meaning |
|------|---------|---------|
| `PROFILE_STAGES` | `True` | Print `[PROFILE]` stage timings |
| `SKIP_QC_PLOTS` | `True` | Skip matplotlib QC PNGs (Excel/TIFs still written) |
| `USE_SLOPE_CACHE` | `True` | Cache AOI-aligned slope stack as `.npz` for re-runs |

## Verify

1. Run original/local baseline and `_opt.py` for the same AOI (e.g. Tanzania).
2. Compare Excel sheets `Top_Soil_Layers` and `Sub_Soil_Layers` (keys exact; weights within float noise).
3. Compare `*_hwsd_smu.tif` and `*_slope_class.tif` pixel values.
4. Compare wall times from the footer: `Total computational time: … seconds`.

## Rollback

Delete or ignore `_opt.py` and keep using `PyAEZ_Module4_Step1_Soil_Slope_v3.py`.

## Not included (on purpose)

- Original v3 (you already have it)
- Step 2 script
- Input data / Tanzania run outputs
- Local baseline/profiling helper scripts

Details: [docs/OPTIMIZATION_NOTES.md](docs/OPTIMIZATION_NOTES.md) · Plain-language: [docs/OPTIMIZATION_NOTES_SIMPLE.md](docs/OPTIMIZATION_NOTES_SIMPLE.md)
