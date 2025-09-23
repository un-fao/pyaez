Modifications
1.  All reduction factors will be externally imported from excel sheets instead of providing
    python scripts.
2.  All reduction factors from excel sheets are recorded as python dictionaries. Algorithm will be the same as
    previous version. But the access of variables will be heavily depending on pandas incorporation and dictionaries.
3.  New soil evaluation method from GAEZ v5 is implemented.
4.  Soil property ratings with zero soil attributes will be set as 100 by default.
5.  Added AWC calculation.

"""

# Script Overview:
# This Python script imports raster and vector spatial datasets along with an Excel file containing soil attributes.
# It processes and extracts relevant soil data within a specified Area of Interest (AOI), cleans and transforms the data,
# and generates an output Excel files for soil layers.
# 2025 (May): Mukaratirwa, RutendoTadiwa (NSLD), Agro-Ecological Zoning Programming Specialist, Rutendo.Mukaratirwa@fao.org
# 2025 (July): Asgharinia, Shahla (NSLD), Crop and Soil Spatial Modelling Specialist, shahla.asgharinia@fao.org

