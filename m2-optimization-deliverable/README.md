# PyAEZ Module 2 — Optimization Package

Faster crop simulation with **identical outputs** to the original `pyaez`. Validated: Malawi rice / 2080 / irrigated / `step_doy=1` — `max_abs_diff = 0`, ~112 s → ~18.5 s at 4 workers (reference PC).

## Package contents

| Path | Purpose |
|------|---------|
| `pyaez/*.py` (4 files) | Drop-in replacements |
| `docs/OPTIMIZATION_NOTES.md` | What changed and why |

## Install

1. Back up `pyaez/` → `pyaez_backup`
2. Copy these four files into your project `pyaez/` (overwrite):
   - `CropSimulation.py`
   - `ThermalScreening.py`
   - `BioMassCalc.py`
   - `CropWatCalc.py`
3. In **Notebook 2** (or your script), add **before** crop simulation:

   ```python
   aez.setParallelJobs(4)   # use 2 on dual-core laptops; 1 for sequential
   ```

4. When using multiple workers, set once per session (notebook or terminal):

   ```python
   import os
   os.environ["OMP_NUM_THREADS"] = "1"
   os.environ["NUMBA_NUM_THREADS"] = "1"
   ```

No other notebook or `data_input` changes are required for the library upgrade.

## Verify

Run **NB2_CropSimulation-user** as usual. The irrigated simulation cell already prints elapsed time (`Simulation took X seconds.`). Outputs (yield maps, fc1/fc2, etc.) should match the pre-upgrade run; only runtime should improve.

To compare speed: note the timing before and after copying the four files (same machine, same `step_doy`, restart kernel).

## Rollback

```powershell
Remove-Item -Recurse pyaez; Copy-Item -Recurse pyaez_backup pyaez
```

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `setParallelJobs` not found | Confirm all four files were copied |
| Slow on 2-core laptop | `aez.setParallelJobs(2)` or `1` |
| High CPU, little speedup | Set `OMP_NUM_THREADS=1` and `NUMBA_NUM_THREADS=1` |

Details: [docs/OPTIMIZATION_NOTES.md](docs/OPTIMIZATION_NOTES.md)
