# Module 4 Step 1 Optimizations — Technical Summary

**Goal:** Same science, faster AOI soil/slope table builds. **Scope:** Step 1 only (`_opt.py` copy).  
**Validated:** Tanzania — Excel sheets `Top_Soil_Layers` / `Sub_Soil_Layers` match baseline (keys exact; weight columns within ~1e-14 float noise). SMU and dominant-slope GeoTIFFs bit-identical.

**Runtime (Tanzania, this machine):** baseline **64.77 s** → optimized **22.37 s** (~**2.9×**).

## Bottleneck (profiled)

| Stage | Baseline | Optimized |
|-------|----------|-----------|
| Load GAUL boundaries | 36.00 s | 0.82 s |
| Align 10 slope-share RSTs | 1.68 s | 2.70 s (cold + write cache) |
| `build_smu_slope_distribution` | 6.36 s | 0.11 s |
| HWSD layers load | 1.33 s | 1.28 s |
| Excel export | (in total) | 14.60 s |
| **Wall time** | **64.77 s** | **22.37 s** |

Largest win: filtered GAUL read. Second: vectorized SMU slope shares. Slope align was already fast with GDAL; windowed reads + NPZ cache help repeated AOI runs more than the first cold pass.

## Principles

- No change to weight formulas, nearest resampling, SMU>0 mask, merge keys, or `CODE` composition
- Do not use `GAEZ-V5.SLOPE-MED.tif` (same as v3)
- Keep original `PyAEZ_Module4_Step1_Soil_Slope_v3.py` unchanged for FAO comparison
- Step 2 left untouched

## Changes in `_opt.py`

| # | Change |
|---|--------|
| 1 | `BASE_DIR` = script parent; `input data/` + `Outputs/optimized/` |
| 2 | Read `HWSD2_LAYERS_2.csv`; filter to AOI SMUs before joins |
| 3 | GAUL `read_file(..., where=...)` for selected AOI only |
| 4 | Windowed AOI-bbox reads of global slope RSTs (nearest unchanged) |
| 5 | Cache AOI-aligned 10-band slope stack as `.npz` (`USE_SLOPE_CACHE`) |
| 6 | Vectorize SMU slope distribution with `np.bincount` |
| 7 | Filter SMUs while parsing `sluslope30_hwsd2v10.dat` |
| 8 | `SKIP_QC_PLOTS=True` (optional matplotlib skip) |
| 9 | Stage timers via `PROFILE_STAGES` |

## How to run

```text
conda activate pyaez
python PyAEZ_Module4_Step1_Soil_Slope_v3_opt.py
```

Flags near the top of the script: `PROFILE_STAGES`, `SKIP_QC_PLOTS`, `USE_SLOPE_CACHE`.

## Not in scope

- Step 2 SQ1–SQ7 / SSR logic
- Modules 1–3
- Changing scientific formulas or sheet schemas Step 2 needs

## One-line summary (for reports)

> Module 4 Step 1 was sped up by filtered GAUL loading, early SMU filtering, windowed slope RST reads with optional cache, and vectorized SMU slope shares—same Tanzania Excel/TIF outputs as v3, ~2.9× faster on this machine.
