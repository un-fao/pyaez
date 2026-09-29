# PyAEZ

[![PyAEZ](https://img.shields.io/badge/PyAEZ-open--source-orange.svg)](#pyaez)
[![Python](https://img.shields.io/badge/Python-3.x-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Developed by](https://img.shields.io/badge/Developed%20by-FAO-116AAB.svg)](https://www.fao.org/geospatial/data-and-tools/tools/python-package-for-agro-ecological-zoning-(pyaez)/en)

PyAEZ is an open and modular Python implementation of the Agro-Ecological
Zoning (AEZ) framework.

It provides transparent, code-based workflows for reproducing, testing,
adapting and extending AEZ assessments. PyAEZ is designed to make AEZ
calculations more accessible, reproducible and adaptable, particularly for
national and subnational applications.

PyAEZ integrates climate, crop, soil and terrain information to simulate
crop production conditions, estimate attainable yields, identify
environmental constraints and assess the agro-ecological suitability of
land for selected crops and production systems.

PyAEZ can use global input datasets or be adapted to national and subnational
data. These may include national or local soil maps and surveys,
meteorological station or gridded climate data, crop calendars and locally
calibrated crop parameters, digital elevation and terrain data, and
user-defined agricultural production and management scenarios.

This flexibility allows users to construct AEZ assessments that reflect
local data availability, agricultural conditions and planning requirements.
Using alternative datasets or modifying calculation procedures requires a
good understanding of the scientific assumptions, data requirements and
sequential logic of the AEZ methodology.

---

## What is Agro-Ecological Zoning?

The Agro-Ecological Zoning (AEZ) methodology builds on the FAO Framework
for Land Evaluation and several decades of methodological development by
the Food and Agriculture Organization of the United Nations (FAO), the
International Institute for Applied Systems Analysis (IIASA) and
collaborating institutions.

AEZ provides a framework for land resource inventory and appraisal and for
assessing the suitability and productivity of land for agricultural uses.

AEZ integrates information on:

- climate and agro-climatic conditions;
- crop characteristics and requirements;
- water availability and crop water requirements;
- soil characteristics and limitations;
- terrain and slope conditions; and
- agricultural production systems and management assumptions.

These data are evaluated through a sequence of crop-specific calculations
to estimate potential and attainable yields and assess the suitability of
land for agricultural production.

---

## PyAEZ framework

The current PyAEZ development framework follows the modular and sequential
approach of the AEZ methodology. The main components are organized from
Module 0 through Module V, together with supporting utilities.

| Module | Name | Main purpose |
|---|---|---|
| Module 0 | Inputs and Scenario Preparation | Prepares climate, crop, soil, terrain and other user-defined input data. |
| Module I | Climate Regime | Calculates agro-climatic indicators, temperature regimes, moisture regimes and growing-period characteristics. |
| Module II | Crop Simulation | Simulates crop calendars, biomass production, crop water requirements and attainable yields. |
| Module III | Climatic Constraints | Applies crop-specific climatic constraints and climate-related yield reductions. |
| Module IV | Soil and Terrain Constraints | Assesses soil quality, soil suitability, slope and terrain effects on crop production. |
| Module V | Agro-Ecological Suitability | Integrates climate, crop, soil and terrain results into final suitability and yield outputs. |
| Utilities | Supporting Tools | Provides preprocessing, validation, raster processing and module-specific utility functions. |

The modules form an interconnected workflow in which outputs from one stage
may provide inputs to subsequent stages. For a complete AEZ assessment, the
main modules should therefore generally be executed sequentially.

The organization and functionality of the modules may evolve as PyAEZ is
further developed.

---

## Main calculation components

PyAEZ includes or supports calculations related to:

- climate regimes and length of growing period;
- reference and crop evapotranspiration;
- crop water balance and crop water requirements;
- crop-calendar simulation;
- biomass and attainable-yield estimation;
- agro-climatic constraints;
- soil quality and crop-specific soil suitability;
- terrain and slope constraints; and
- integrated agro-ecological suitability.

Several calculation components can also be used independently without
running the complete PyAEZ workflow.

---

## Using PyAEZ

PyAEZ is implemented in Python and can be used through Python scripts and
Jupyter Notebooks.

Jupyter-based workflows provide an interactive environment in which users
can examine input data, execute individual calculation steps, inspect
intermediate results and adapt parameters for specific applications.

The use of Jupyter Notebooks also allows PyAEZ workflows to be implemented
in cloud-based notebook environments such as Google Colab, supporting
training and capacity-development activities where local software
installation may be impractical.

The current development workflow includes notebooks corresponding to the
main stages of the assessment:

```text
M0_Inputs.ipynb
M1_ClimateRegime.ipynb
M2_CropSimulation.ipynb
M3_ClimaticConstraints.ipynb
M4_SoilTerrainConstraints.ipynb
M5_AgroecologicalSuitability.ipynb
```

The exact notebooks, interfaces and module organization may change between
PyAEZ releases. Users should refer to the documentation accompanying the
release being used.

---

## Repository structure

The repository is organized around the PyAEZ source code, calculation
modules, utilities, documentation and example workflows.

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
│   │
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

This structure reflects the current development organization and may change
in future releases.

---

## Scope and limitations

PyAEZ is a scientific modelling framework. Results depend on the quality,
representativeness and spatial and temporal resolution of the input
datasets, crop parameters, production assumptions and methodological
configuration used in an assessment.

Users should consider:

- the quality and representativeness of climate, soil and terrain data;
- the appropriateness of crop parameters for local conditions;
- assumptions regarding water supply and agricultural management;
- spatial and temporal resolution;
- calibration and validation against observations where possible; and
- the applicability of individual algorithms to the environmental
  conditions being assessed.

PyAEZ outputs should be interpreted together with the assumptions and
limitations of the underlying datasets and methodologies.

---

## Citation

When using PyAEZ in a publication, report or technical assessment, please
cite the software and identify the version used.

For the recommended software citation, see the [`CITATION.cff`](CITATION.cff)
file or use the **Cite this repository** option provided by GitHub.

---

## Development history and acknowledgements

PyAEZ originated from efforts to make Agro-Ecological Zoning (AEZ)
methodologies more accessible for national and subnational applications
through an open and reproducible Python-based modelling environment.

PyAEZ has benefited from collaboration between FAO, the Asian Institute of
Technology - Geo-informatics Center (AIT-GIC), the International Institute
for Applied Systems Analysis (IIASA), and Mississippi State University (MSU).

The current development of PyAEZ at FAO has benefited from contributions
from Shahla Asgharinia, Dario Spiller, Rutendo Tadiwa Mukaratirwa,
Federica Chiozza, Livia Peiser, Gianluca Franceschini, Matieu Henry,
and Beau Damen.

PyAEZ has also benefited from technical contributions and advice from
Günther Fischer (IIASA) and Freddy Nachtergaele.

The current development also includes the integration of PyAEZ and CAVApy,
in collaboration with the FAO OCB Division.

Further information on the earlier implementation and its contributors is
available in the [AIT-GIC PyAEZ repository](https://github.com/gicait/PyAEZ).

FAO acknowledges all institutions, projects and individuals that have
contributed to the development and evolution of PyAEZ.

The current development of PyAEZ at FAO is financially supported by the
SoilFER programme.

---

## License

See the repository [`LICENSE`](LICENSE) file for the applicable terms
governing the use, modification and distribution of PyAEZ.

---

## Major AEZ references

1. de Wit, C.T. 1965. *Photosynthesis of Leaf Canopies*. Agricultural
   Research Report No. 663. PUDOC, Wageningen.

2. FAO. 1992. *CROPWAT: A Computer Program for Irrigation Planning and
   Management*. FAO Irrigation and Drainage Paper No. 46. Rome.

3. FAO. 1998. *Crop Evapotranspiration: Guidelines for Computing Crop Water
   Requirements*. FAO Irrigation and Drainage Paper No. 56. Rome.

4. FAO. 2017. *Final Report: National Agro-Economic Zoning for Major Crops
   in Thailand*.

5. Fischer, G., van Velthuizen, H., Shah, M. and Nachtergaele, F. 2002.
   *Global Agro-Ecological Assessment for Agriculture in the 21st Century:
   Methodology and Results*. IIASA Research Report RR-02-02. Laxenburg,
   Austria.

6. Fischer, G., Nachtergaele, F.O., van Velthuizen, H., Chiozza, F.,
   Franceschini, G., Henry, M., Muchoney, D. and Tramberend, S. 2021.
   *Global Agro-Ecological Zones v4 - Model Documentation*. FAO, Rome.

7. Monteith, J.L. 1965. Evapotranspiration and the environment. In:
   *The State and Movement of Water in Living Organisms*, pp. 205-234.

8. Monteith, J.L. 1981. Evapotranspiration and surface temperature.
   *Quarterly Journal of the Royal Meteorological Society*, 107: 1-27.

9. FAO and IIASA. 2026. *Soil Suitability Assessment Procedures for
   Conventional and Organic Farming Systems Used in Global
   Agro-Ecological Zoning Version 5*. Rome and Laxenburg, Austria.

---

## 🔗 Useful Links

- 🐍 **FAO PyAEZ Website:** [Python package for Agro-Ecological Zoning (PyAEZ)](https://www.fao.org/geospatial/data-and-tools/tools/python-package-for-agro-ecological-zoning-(pyaez)/en)
- 📖 **GAEZ v5 Documentation:** [GAEZ v5 Wiki](https://github.com/un-fao/gaezv5/wiki)
