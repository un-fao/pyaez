# Module 4 Step 1 — What I Optimized (Simple Summary)

## In one sentence

I made **Step 1** build the Tanzania soil/slope Excel and maps **about 3× faster**, without changing the scientific results. **Step 2 was not changed.**

## What Step 1 does (plain English)

Step 1 gathers soil and slope information for a country (or other area) and writes:

- An Excel workbook (soil layers + slope weights)
- Two map files (soil units and dominant slope class)

Step 2 later uses those files. I only sped up Step 1.

## What stayed the same

- Same formulas for weights and slope classes  
- Same Excel sheet layout Step 2 expects  
- Same map values (checked against a baseline Tanzania run)  
- The original FAO script (`…_v3.py`) is unchanged so FAO can compare  

## What I made faster

Most of the old runtime was spent loading and handling **much larger data than needed** for one country.

1. **Country boundaries** — read only Tanzania from the world boundary file, instead of loading the whole world first  
2. **Soil attribute table** — keep only the soil units that appear in Tanzania  
3. **Slope maps** — work on the Tanzania area of the global slope files, not the entire globe each time; optionally save that result so a second run is quicker  
4. **Summaries by soil unit** — replace a slow loop with a faster calculation that gives the same numbers  
5. **Optional skip of QC pictures** — Excel and maps still write; the preview PNGs can be turned off  

## Speed (my test PC)

| Version | Time (Tanzania) |
|---------|-----------------|
| Before (baseline) | about **65 seconds** |
| After (optimized) | about **22 seconds** |

Roughly **2.9× faster**. Times will vary by computer and country size.

## How to use it

1. Keep the original FAO Step 1 script.  
2. Add `PyAEZ_Module4_Step1_Soil_Slope_v3_opt.py` next to it.  
3. Run the `_opt` script with your usual Module 4 inputs.  

More detail (for developers): [OPTIMIZATION_NOTES.md](OPTIMIZATION_NOTES.md).

## What this package is *not*

- Not a rewrite of Module 4 science  
- Not an update to Step 2 (waiting on FAO’s final Step 2)  
- Not a change to Modules 1–3  
