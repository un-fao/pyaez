[![PyAEZ](https://img.shields.io/badge/PyAEZ-3.0%20development-orange.svg)](#pyaez-version-30)
[![Release status](https://img.shields.io/badge/release-forthcoming-yellow.svg)](#release-status)
[![Python](https://img.shields.io/badge/Python-3.14-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Developed by](https://img.shields.io/badge/Developed%20by-FAO%20NSL-116AAB.svg)](https://www.fao.org/land-water/en/)

## Modules

The current PyAEZ Version 3.0 repository contains Module 0 through Module V.

| Module | Name | Main purpose |
|---|---|---|
| Module 0 | Inputs and Scenario Preparation | Prepares climate, crop, soil, terrain and other user-defined input data. |
| Module I | Climate Regime | Calculates agro-climatic indicators, temperature regimes, moisture regimes and growing-period characteristics. |
| Module II | Crop Simulation | Simulates crop calendars, biomass production, crop water requirements and attainable yields. |
| Module III | Climatic Constraints | Applies crop-specific climatic constraints and climate-related yield reductions. |
| Module IV | Soil and Terrain Constraints | Assesses soil quality, soil suitability, slope and terrain effects on crop production. |
| Module V | Agro-Ecological Suitability | Integrates climate, crop, soil and terrain results into final suitability and yield outputs. |
| Utilities | Supporting Tools | Provides preprocessing, validation, raster processing and module-specific utility functions. |

The broader AEZ framework may also include economic suitability analysis.
This component is not currently included in the PyAEZ Version 3.0 workflow.

---

## Step-by-Step Workflow

The following Jupyter notebooks provide worked examples for the main
modules:

1. `M0_Inputs.ipynb`
2. `M1_ClimateRegime.ipynb`
3. `M2_CropSimulation.ipynb`
4. `M3_ClimaticConstraints.ipynb`
5. `M4_SoilTerrainConstraints.ipynb`
6. `M5_AgroecologicalSuitability.ipynb`

The notebooks should generally be run sequentially because outputs from
one module may be required by the following modules.

---

## Repository Structure

```text
Python-AEZ-FAO/
│
├── PyAEZ3.0/
│   ├── module0/
│   ├── module1/
│   ├── module2/
│   ├── module3/
│   ├── module4/
│   ├── module5/
│   │
│   ├── data_input/
│   ├── data_output/
│   │
│   ├── utilities/
│   │   ├── module0/
│   │   ├── module1/
│   │   ├── module2/
│   │   ├── module3/
│   │   ├── module4/
│   │   └── module5/
│   ├── pyaez/
│   ├── docs/
│   ├── website_docs/
│   └── README.md
│
├── LICENSE
├── CITATION.cff
├── environment_pyaez.yml
└── README.md
```

## Citation

When using PyAEZ Version 3.0 in a publication, report or technical
assessment, cite the repository and specify the release tag or Git commit
used.

### Recommended citation




## License
See the repository LICENSE file for the applicable terms.
