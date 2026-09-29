[![PyAEZ](https://img.shields.io/badge/PyAEZ-3.0%20development-orange.svg)](#pyaez-version-30)
[![Release status](https://img.shields.io/badge/release-forthcoming-yellow.svg)](#release-status)
[![Python](https://img.shields.io/badge/Python-3.14-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Developed by](https://img.shields.io/badge/Developed%20by-FAO%20NSL-116AAB.svg)](https://www.fao.org/land-water/en/)

# PyAEZ Version 3.0

PyAEZ is an open and modular Python implementation of the Agro-Ecological
Zoning (AEZ) framework.

It provides transparent, code-based workflows for reproducing, testing,
adapting and extending AEZ assessments. Each module can be inspected,
modified, debugged and independently validated.

PyAEZ can use global datasets or be adapted to national and subnational
information, including digital soil maps, local soil surveys, updated
climate datasets, crop parameters and user-defined production scenarios.

Applying alternative datasets or modifying the calculation procedures
requires a good understanding of the scientific assumptions, input
requirements and sequential logic of the AEZ methodology.

---

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

## Main Calculation Components

The framework includes or supports the following main calculations:

- climate-regime analysis;
- length-of-growing-period assessment;
- reference and crop evapotranspiration;
- crop water-balance calculations;
- biomass and attainable-yield estimation;
- agro-climatic constraint assessment;
- soil-quality and soil-suitability rating;
- terrain and slope-impact assessment; and
- final agro-ecological suitability classification.

Several calculation components can also be used independently without
running the complete PyAEZ workflow.

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

## Funding

## Major AEZ References
1. de Wit, C.T. 1965. Photosynthesis of Leaf Canopies. Agricultural
Research Report No. 663. PUDOC, Wageningen.
2. FAO. 1992. CROPWAT: A Computer Program for Irrigation Planning and
Management. FAO Irrigation and Drainage Paper No. 46. Rome.
3. FAO. 1998. Crop Evapotranspiration: Guidelines for Computing Crop Water
Requirements. FAO Irrigation and Drainage Paper No. 56. Rome.
4. FAO. 2017. Final Report: National Agro-Economic Zoning for Major Crops
in Thailand.
5. Fischer, G., van Velthuizen, H., Shah, M. and Nachtergaele, F. 2002.
Global Agro-Ecological Assessment for Agriculture in the 21st Century:
Methodology and Results. IIASA Research Report RR-02-02.
6. Monteith, J.L. 1965. Evapotranspiration and the environment. In:
The State and Movement of Water in Living Organisms, pp. 205–234.
7. Monteith, J.L. 1981. Evapotranspiration and surface temperature.
Quarterly Journal of the Royal Meteorological Society, 107: 1–27.
8. FAO and IIASA. 2026. Soil Suitability Assessment Procedures for
Conventional and Organic Farming Systems Used in Global
Agro-Ecological Zoning Version 5. Rome and Laxenburg, Austria.
https://doi.org/10.4060/cd8206en
