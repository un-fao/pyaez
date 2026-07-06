# Module 2 Optimizations — Technical Summary

**Goal:** Same science, faster runs. **Scope:** four `pyaez` files only.  
**Validated:** Malawi rice, 2080, irrigated, `step_doy=1` — all six output layers match baseline (`max_abs_diff = 0`).

## Bottleneck

Per valid pixel: up to **365 planting dates** × (thermal screening + biomass + water balance). Baseline repeated heavy setup on every cycle and used one CPU core.

## Principles

- No equation or threshold changes
- Hoist work that is identical across all 365 cycles (once per pixel or once per run)
- Parallelize **across pixels** (process pool), not inside water-balance Numba loops
- Backward compatible: `setParallelJobs(1)` ≈ old sequential behavior

## Changes by file

| # | File | Change |
|---|------|--------|
| 9 | `CropSimulation.py` | `setParallelJobs(n)`, `ProcessPoolExecutor` over valid pixels |
| 6 | `CropSimulation.py` | Duplicate 1y→2y climate once per run, not per pixel |
| 4 | `CropSimulation.py` | Vectorized valid-pixel mask; cleaner pixel driver |
| 2 | `CropSimulation.py` | Pre-sized cycle arrays instead of `np.append` |
| 1 | `CropSimulation.py` | Skip crop-rule profiling when `setCropSpecificRule=False` |
| 3 | `ThermalScreening.py` | Cache year-level crop-rule context; per-cycle window only |
| 5 | `BioMassCalc.py` | Cache latitude splines + 731-day radiation LUT |
| 8 | `BioMassCalc.py` | Hoist LAI → growth-rate multiplier (constant per pixel) |
| 7 | `CropWatCalc.py` | Hoist stage day bounds + rooting depth (Numba helpers) |

**Largest wins:** #9 (multi-core) and #3/#5/#7 (stop redoing work 365× per pixel).

## Not in scope

- No vectorized 365-day loop; no Module 0/1 changes
- Existing Numba kernels (`getReductionFactorNumba`, `calculateMoistureLimitedYieldNumba`, etc.) reused, not rewritten
- Numba `parallel=True` on water balance was **not** used (known to slow that kernel)

## API

```python
aez.setParallelJobs(4)  # default workers
aez.simulateIrrigatedCropCycle(..., n_jobs=4)  # optional override; -1 = all cores
```

When `n_jobs > 1`: set `OMP_NUM_THREADS=1` and `NUMBA_NUM_THREADS=1`.

## Tuning

| Hardware | `n_jobs` |
|----------|----------|
| 1 core | 1 (still faster than baseline from opts #2–#8) |
| 2 cores | 2 |
| 4+ cores | 4 |

## Known limits

- `simulateIrrigatedCropSpecific()` — pre-existing bug, unchanged
- Large countries / low RAM — reduce workers if memory is tight

## One-line summary (for reports)

> Module 2 was sped up by caching per-pixel invariants, hoisting climate duplication, and process-parallel pixel simulation—four files only, identical Malawi rice outputs, ~6× faster at four workers.
