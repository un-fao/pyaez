""""
PyAEZ version 2.3(Apr 2024)
2022- 2023: Swun Wunna Htet, Kittiphon Boonma
2023 (Dec) : Swun Wunna Htet
2024 (Apr) : Swun Wunna Htet, Dwijendra Das

Modification:

1.  Based of GAEZ appendices, the two different lookup tables of reduction factors
    depending on the annual mean temperatures of >= 20 and < 10 deg C are added.
    With this new tables, new fc3 factors are calculated separately for rainfed and
    irrigated conditions.
    The algorithm will check the annual mean temperature and assess the value from 
    the respective look-up fc3 factor to apply to yield.
2. Added missing logic of linear interpolation for pixels with annual mean temperature between
    10 and 20 deg Celsius to extract fc3 constraint factor.
3. Adding missing logic of linear interpolation for wetness-day-specific agro-climatic constraints
4. Excel sheets of agro-climatic constraint factors are required to provide into the system instead of python file.
5. A new object class (ProcessExcelWrapper) is added for excel sheet creation for agro-climatic constraints.
    
"""
import numpy as np
import pandas as pd
from pyaez.ETOCalc import calculateETONumba
from pyaez.UtilitiesCalc import generateLatitudeMap, interpMonthlyToDaily, averageDailyToMonthly
import warnings
import os
import matplotlib.pyplot as plt

try:
    from osgeo import gdal
except:
    import gdal
    
warnings.filterwarnings('ignore')

np.round_ = np.round

class ClimaticConstraints(object):

    def __init__(self, lat_min, lat_max, elevation, mask = None, no_mask_value = None):
        """Calling object class of Climate Constraints. Providing minimum and maximum latitudes, and mask layer.
        
        Args:
            lat_min (float): Minimum latitude [Unit: Decimal Degrees]
            lat_max (float) : Maximum latitude [Unit: Decimal Degrees]
            elevation (float/integer): elevation [Unit: meters]
            mask [integer]: mask layers [binary, 0/1]
        Return:
            None.
        """
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.elevation = elevation
        self.set_daily = False
        self.set_monthly = False
        
        self.im_height = elevation.shape[0]
        self.im_width = elevation.shape[1]
        self.latitude = generateLatitudeMap(lat_min, lat_max, self.im_height, self.im_width)
        self.leap_year = False
        self.set_mask = False

        if mask is not None:
            self.im_mask = mask
            self.no_mask_value = no_mask_value
            self.set_mask = True

    def setClimateData(self, min_temp, max_temp, wind_speed, short_rad, rel_humidity, precip):
        """Load the DAILY or MONTHLY climatic data into Module III object class.
        
        Args:
            min_temp (3D NumPy, float): Minimum temperature [Unit: Celsius]
            max_temp (3D NumPy, float): Maximum temperature [Unit: Celsius]
            wind_speed (3D NumPy, float): Windspeed at 2m altitude [Unit: m/s]
            short_rad (3D NumPy, float): Radiation [Unit: W/m2]
            rel_humidity (3D Numpy, float): Relative humidity [Unit: decimal percentage]
            precipitation (3D Numpy, float): Precipitation [Unit: mm/day]
        Return:
            None.
        """
        rel_humidity[rel_humidity > 0.99] = 0.99
        rel_humidity[rel_humidity < 0.05] = 0.05
        short_rad[short_rad < 0] = 0
        wind_speed[wind_speed < 0] = 0

        # Time dimension checkpoint
        doy = None
        if np.all(min_temp.shape[2] ==12 and max_temp.shape[2] ==12 and wind_speed.shape[2] ==12
                and short_rad.shape[2] ==12 and rel_humidity.shape[2] ==12 and precip.shape[2] ==12):
            self.set_monthly = True
            doy = 365
        elif np.all(min_temp.shape[2] ==365 and max_temp.shape[2] ==365 and wind_speed.shape[2] ==365
                and short_rad.shape[2] ==365 and rel_humidity.shape[2] ==365 and precip.shape[2] ==365):
            doy = 365
        elif np.all(min_temp.shape[2] ==366 and max_temp.shape[2] ==366 and wind_speed.shape[2] ==366
                and short_rad.shape[2] ==366 and rel_humidity.shape[2] ==366 and precip.shape[2] ==366):
            doy = 366
            self.leap_year = True
        else:
            raise ValueError('Time Dimension of climate data must be 12, 365 or 366. Please check your input data.')

        self.meanT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.totalPrec_daily = np.zeros((self.im_height, self.im_width, doy))
        self.minT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.maxT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.shortrad_daily = np.zeros((self.im_height, self.im_width, doy))
        self.wind_daily = np.zeros((self.im_height, self.im_width, doy))
        self.rel_humidity_daily = np.zeros((self.im_height, self.im_width, doy))
        self.eto_daily = np.zeros((self.im_height, self.im_width, doy))
        self.shortrad_daily_MJm2day = np.zeros((self.im_height, self.im_width, doy))

        self.min_T = np.zeros((self.im_height, self.im_width))
        self.months_P_gte_eto = np.zeros((self.im_height, self.im_width))
        monthly_precip = np.zeros((self.im_height, self.im_width, 12))
        monthly_eto = np.zeros((self.im_height, self.im_width, 12))

        if self.set_monthly:

            for i_row in range(self.im_height):
                for i_col in range(self.im_width):

                    if self.set_mask:
                        if self.im_mask[i_row, i_col] == self.no_mask_value:
                            continue
                    
                    Tm = (max_temp[i_row, i_col, :] + min_temp[i_row, i_col, :])/2
                    self.meanT_daily[i_row, i_col, :] = interpMonthlyToDaily(Tm, 1, doy)
                    self.totalPrec_daily[i_row, i_col, :] = interpMonthlyToDaily(precip[i_row, i_col,:], 1, doy, no_minus_values=True)
                    self.minT_daily[i_row, i_col, :] = interpMonthlyToDaily(min_temp[i_row, i_col,:], 1, doy)
                    self.maxT_daily[i_row, i_col, :] = interpMonthlyToDaily(max_temp[i_row, i_col,:], 1, doy)
                    self.shortrad_daily[i_row, i_col, :] = interpMonthlyToDaily(short_rad[i_row, i_col,:], 1, doy, no_minus_values=True)
                    self.wind_daily[i_row, i_col, :] = interpMonthlyToDaily(wind_speed[i_row, i_col,:], 1, doy, no_minus_values=True)
                    self.rel_humidity_daily[i_row, i_col, :] = interpMonthlyToDaily(rel_humidity[i_row, i_col,:], 1, doy, no_minus_values=True)

                    monthly_precip[i_row, i_col, : ] = precip[i_row, i_col, :]
                    self.min_T[i_row, i_col] = np.nanmin(Tm)
                
        else:
            self.meanT_daily = (min_temp+ max_temp)/2
            self.totalPrec_daily = precip
            self.minT_daily = min_temp
            self.maxT_daily = max_temp
            self.shortrad_daily = short_rad
            self.wind_daily = wind_speed
            self.rel_humidity_daily = rel_humidity

            for i_row in range(self.im_height):
                for i_col in range(self.im_width):
                    monthly_precip[i_row, i_col, : ]=  averageDailyToMonthly(self.totalPrec_daily[i_row, i_col, :], self.leap_year)        
                    monthly_Tm = averageDailyToMonthly(self.meanT_daily[i_row, i_col, :], self.leap_year)
                    self.min_T[i_row, i_col] = np.nanmin(monthly_Tm)
        # ET0 calculation
        self.shortrad_daily_MJm2day = (self.shortrad_daily*3600*24)/1000000 # convert w/m2 to MJ/m2/day
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):
                self.eto_daily[i_row, i_col, :] = calculateETONumba(1, doy, self.latitude[i_row, i_col], self.elevation[i_row, i_col], 
                                                                    self.minT_daily[i_row, i_col, :], self.maxT_daily[i_row, i_col, :], 
                                                                    self.wind_daily[i_row, i_col, :], self.shortrad_daily_MJm2day[i_row, i_col, :],
                                                                      self.rel_humidity_daily[i_row, i_col, :])
                
                monthly_eto[i_row,i_col,:] = averageDailyToMonthly(self.eto_daily[i_row,i_col,:], self.leap_year)
        
        # counting months with monthly precipitation >= monthly ET0
        self.months_P_gte_eto = np.sum(monthly_precip >= monthly_eto, axis = 2)
    

    def setReductionFactors(self, file_path):
        """ Load the agro-climatic reduction factors for either rainfed or irrigated conditions.

        Args:
            file_path (String): The directory file path of excel sheet in xlsx format storing agro-climatic reduction factor.
                                The excel must contain three sheets namely: mean>20, mean<10 and lgpt10.
        Return: 
            None.
        """
        main = pd.read_excel(file_path, sheet_name=None)

        if main['lgpt10'].isnull().values.any()==True or main['mean>20'].isnull().values.any()==True or  main['mean<10'].isnull().values.any()==True:
            print('Missing values of reduction factor detected. Excel sheets with no null-values required')
            del(main)
        
        else:
            self.gte20 = main['mean>20']
            self.lt10 = main['mean<10']
            self.lgpt10 = main['lgpt10']
            del(main)
            
    def setReductionFactors_fromM0(self, condition_type):
        """
        Load agro‑climatic reduction factors previously saved in HDF5 format.
    
        This method reads the reduction factor tables for the three agro‑climatic
        constraint categories (gte20, lt10, lgpt10) from an HDF5 file created by the
        AEZ pre‑processing pipeline. Each table is stored as a pandas DataFrame.
    
        Parameters
        ----------
        file_path : str
            Path to the HDF5 (.h5) file containing the stored reduction factor
            DataFrames. The file must contain the following keys:
            - "gte20"   : reduction factors for mean ≥ 20% probability
            - "lt10"    : reduction factors for mean < 10% probability
            - "lgpt10"  : reduction factors for length of growing period > 10 months
    
        Notes
        -----
        - This function does not perform validation on the structure of the loaded
          DataFrames; it assumes they were saved using the AEZ utility functions.
        - The loaded tables are assigned to the instance attributes:
          `self.gte20`, `self.lt10`, and `self.lgpt10`.
    
        Returns
        -------
        None
            The method updates the object's internal state but does not return data.
    
        Raises
        ------
        KeyError
            If any of the required keys ("gte20", "lt10", "lgpt10") are missing
            from the HDF5 file.
        IOError
            If the file cannot be opened or read.
        """

        file_path = f"./data_input/input_module3/M3_constraints_{self.crop_name}_{self.input_level}_{condition_type}.h5"
        
        # Reading and storing the saved file
        with pd.HDFStore(file_path) as store:
            self.gte20  = store["gte20"]
            self.lt10   = store["lt10"]
            self.lgpt10 = store["lgpt10"]    

    def calculateLGPagc(self, lgp, lgp_equv):
        """ Calculation of adjustted LGP for agro-climatic constraints.
        
        Args:
            lgp (Numerical): Length of Growing Period [Unit: Days]
            lgp_equv (Numerical): Equivalent Length of Growing Periods [Unit: Days]

        Return:
            lgp_agc (int): Adjusted LGP for agro-climatic constraints [Unit: Days]. 
        """

        # Wetness indicator calculation logic referred to GAEZ v4 Model Documentation Pg. 72

        if lgp <= 120:
            lgp_agc= min(120, max(lgp, lgp_equv)) # correct
        
        elif lgp in range(121,210+1):
            lgp_agc = lgp # correct
            
        elif lgp > 210:
            lgp_agc = max(210, min(lgp, lgp_equv)) # correct
        
        else:
            raise ValueError('Something is wrong.')
        
        return lgp_agc
    
    def applyClimaticConstraints(self, yield_input, lgp, lgp_equv, lgpt10, omit_yld_0= False):

        """
        Args:
            yield_input (2D-NumPy Array): Yield map to apply agro-climatic constraint factor. [Unit: kg/ha]
            lgp (2D NumPy Array): Length of Growing Period [Unit: Days]
            lgp_equv (2D NumPy Array): Equivalent Length of Growing Periods [Unit: Days]
            lgpt10 (2D NumPy Array): Thermal Growing Periods at 10 degrees [Unit: Days]
            omit_yld_0 (Boolean): Any zero yield areas will not be calculated. Default is False.
        Return:
            None.
        """

        self.adj_yield = np.zeros((self.im_height, self.im_width), dtype = int)
        original_yld = np.copy(yield_input)
        self.lgp_agc = np.zeros((self.im_height, self.im_width), dtype = int)
        self.fc3 = np.zeros((self.im_height, self.im_width), dtype = np.float16)

        # Middle day of year for each agro-climatic constraints (used for linear interpolation purposes)
        mid_doy = np.array([0, 15,  45,  75, 105, 135, 165, 195, 225, 255, 285, 315, 345, 365]) # total 14 interval points

        for i in range(self.im_height):
            for j in range(self.im_width):

                if self.set_mask:
                    if self.im_mask[i,j] == self.no_mask_value:
                        continue
                
                if omit_yld_0:
                    if original_yld[i,j] == 0:
                        continue

                # for LPG having 365 or 366, either 365+ or 365- will be selected
                self.lgp_agc[i,j] = self.calculateLGPagc(lgp[i,j], lgp_equv[i,j])
                
                if self.lgp_agc[i,j] >=365 and self.months_P_gte_eto[i,j] == 12:
                    gte20 = (self.gte20.drop(columns = ['365-', 'type'])).to_numpy()
                    lt10 = (self.lt10.drop(columns=['365-', 'type'])).to_numpy()

                else:
                    gte20 = (self.gte20.drop(columns = ['365+', 'type'])).to_numpy()
                    lt10 = (self.lt10.drop(columns=['365+', 'type'])).to_numpy()
                

                # Appending zero reduction factor for zero LGPagc
                # gte20 = np.append(0,gte20)
                # lt10 = np.append(0,lt10)
                
                # Annual mean temperature will select the relevant look-up table
                
                # Case I: ann_mean >= 20
                if self.min_T[i,j] >= 20:
                    B_row = 1 - (np.append(0, gte20[0,:]) / 100)
                    C_row = 1 - (np.append(0, gte20[1,:])/100)
                    D_row = 1 - (np.append(0, gte20[2,:])/100)
                
                # Case II: ann_mean <= 10:
                elif self.min_T[i,j] <= 10:
                    B_row = 1 - (np.append(0, lt10[0,:])/100)
                    C_row = 1 - (np.append(0, lt10[1,:])/100)
                    D_row = 1 - (np.append(0, lt10[2,:])/100)
                
                # Case III: ann_mean between 10 and 20. Linear interpolation is applied.
                else:

                    # 'B' constraint row interpolation
                    B_row_10 = np.append(0, lt10[0,:])
                    B_row_20 = np.append(0, gte20[0,:])
                    B_row = np.zeros(B_row_10.shape[0])

                    for e in range(B_row_10.shape[0]):
                        B_row[e] = 1 - ((np.interp(self.min_T[i,j], [10,20], [B_row_10[e], B_row_20[e]]))/ 100)
                    
                    
                    # 'C' constraint row interpolation
                    C_row_10 = np.append(0, lt10[1,:])
                    C_row_20 = np.append(0, gte20[1,:])
                    C_row = np.zeros(B_row_10.shape[0])

                    for e in range(C_row_10.shape[0]):
                        C_row[e] = (1 - (np.interp(self.min_T[i,j], [10,20], [C_row_10[e], C_row_20[e]]))/100)
                    

                    # 'D' constraint row interpolation
                    D_row_10 = np.append(0, lt10[2,:])
                    D_row_20 = np.append(0, gte20[2,:])
                    D_row = np.zeros(B_row_10.shape[0])

                    for e in range(D_row_10.shape[0]):
                        D_row[e] = 1 - ((np.interp(self.min_T[i,j], [10,20], [D_row_10[e], D_row_20[e]]))/100)
                
                
                # 'E' constraint row interpolation
                E_row = np.append(0, self.lgpt10.drop(columns= 'type').iloc[0].to_numpy())
                E_row = 1 - (E_row/100)

                # Start calculation of agro-climatic constraints
                # 1: find agro-climatic factors of its corresponding interval of wetness days 
                # 2: select the most limiting factor amongst 'b', 'c', 'd' and 'e' constraints

                B = np.interp(self.lgp_agc[i,j], mid_doy, B_row)
                C = np.interp(self.lgp_agc[i,j], mid_doy, C_row)
                D = np.interp(self.lgp_agc[i,j], mid_doy, D_row)
                E = np.interp(lgpt10[i,j], mid_doy, E_row)

                self.fc3[i,j] = np.round(np.min([B*C*D, E]), 2)

                self.adj_yield[i,j] = int(np.round(original_yld[i,j] * self.fc3[i,j], 0))


    def getClimateAdjustedYield(self):
        """
        Generate yield map adjusted with agro-climatic constraints.

        Args:
            None.
        Return:
            clim_adj_yld (2D-Numpy Array): Agro-climatic constraint applied yield [Unit: kg/ha].
        """
        return self.adj_yield
    
    def getClimateReductionFactor(self):
        """
        Generates agro-climatic constraint map (fc3) applied to unconstrainted 
        yield.

        Args:
            None.
        Return:
            fc3 (2D-Numpy Array): agro-climatic constraint factor (fc3).
        """
        return self.fc3

    def compute_yield(
            self,
            condition_type: str,
            plot_results: bool = False
        ):
        """
        Compute climate-adjusted yield and reduction factor for rainfed or irrigated
        conditions, selecting the correct parameter file based on the crop name and
        management level (HI/LI).
    
        This function:
          - Loads the baseline yield raster (NB2) for the requested condition,
          - Loads agro-climatic indicators (LGP, LGPT10, LGPEquivalent) from NB1,
          - Selects the appropriate parameter workbook (Module 3) by `crop_name` and `management`,
          - Applies climatic reduction factors and constraints via `clim_con`,
          - Returns the climate-adjusted yield and the overall reduction factor (Fc3),
          - Optionally plots original yield, constrained yield, and Fc3.
    
        Parameters
        ----------
        self : pyaez.ClimaticConstraints.ClimaticConstraints
            Configured climatic constraints object. Must provide:
            `crop_name` attribute and methods:
            `setReductionFactors(file_path)`,
            `applyClimaticConstraints(...)`,
            `getClimateAdjustedYield()`,
            `getClimateReductionFactor()`.
        condition_type : {'rainfed', 'irrigated'}
            Production condition to compute. Determines which baseline yield raster
            is loaded (NB2: `yld_rain` vs `yld_irr`) and which parameter file is used.
        country_name : str
            Country/study-area identifier used in NB1/NB2 filenames.
        year : int or str
            Year used in NB1/NB2 filenames. Accepts either `int` (e.g., 1990) or `str` (e.g., "1990").
        management : {'HI','LI'}, optional
            Parameter set selection:
            - 'HI' → High-input workbooks (Module 3),
            - 'LI' → Low-input workbooks (Module 3).
            Raises `NameError` if not one of {'HI','LI'}.
        plot_results : bool, optional
            If True, displays three plots (original yield, climate-adjusted yield, Fc3).
    
        Returns
        -------
        clim_yield : numpy.ndarray
            Climate-constrained yield map (same shape as input yield).
        fc3 : numpy.ndarray
            Overall climatic reduction factor (Fc3), values typically in [0, 1].
    
        Notes
        -----
        - **File dependencies**:
            * NB2 baseline yield rasters in `data_output/NB2`:
              `{country_name}_{crop_name}_yld_rain_{year}.tif` or `{country_name}_{crop_name}_yld_irr_{year}.tif`
            * NB1 agro-climatic indicators in `data_output/NB1`:
              `{country_name}_LGP_{year}.tif`,
              `{country_name}_LGPt10_{year}.tif`,
              `{country_name}_LGPEquivalent_{year}.tif`
            * Module 3 parameter workbooks in `data_input/input_module3/`:
              crop-specific `.xlsx` files (HI or LI set)
        - **Crop name normalization**: The function derives `crop_key` from `self.crop_name`
          by dropping trailing `_H` or `_L` and lowercasing (e.g., `"Maize_HI"` → `"maize"`).
        - **Array shapes**: All rasters are expected to be co-registered (same extent/resolution/CRS).
        - **Plotting**: Uses Matplotlib; `vmax` is set to the maximum of original vs. constrained yield.
    
        Raises
        ------
        ValueError
            If `condition_type` is not one of {'rainfed','irrigated'} or if `crop_key` is unsupported.
        NameError
            If `management` is neither 'HI' nor 'LI'.
        FileNotFoundError / OSError
            If any of the required rasters or parameter files cannot be opened.
    
        Examples
        --------
        >>> clim_yield, fc3 = compute_yield(
        ...     self,
        ...     condition_type='rainfed',
        ...     country_name='Ghana',
        ...     year=1990,
        ...     management='HI',
        ...     plot_results=True
        ... )
        >>> clim_yield.shape, fc3.min(), fc3.max()
           ((rows, cols), 0.0, 1.0)
        """
        
        crop_name = self.crop_name
        country_name = self.country_name
        year = self.year
        management = self.input_level
        # Set working directory
        work_dir=os.getcwd()
    
        # Validate condition type
        if condition_type not in ["rainfed", "irrigated"]:
            raise ValueError("condition_type must be 'rainfed' or 'irrigated'")
    
        # Apply reduction factors and constraints
        self.setReductionFactors_fromM0(condition_type)

        # Load yield map
        if condition_type == "rainfed":
            yield_map = gdal.Open(os.path.join(work_dir,f'data_output/NB2/{country_name}_{crop_name}_yld_rain_{year}.tif')).ReadAsArray()
        else:
            yield_map = gdal.Open(os.path.join(work_dir,f'data_output/NB2/{country_name}_{crop_name}_yld_irr_{year}.tif')).ReadAsArray()

        # Load agro-climatic indicators
        lgp = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGP_{year}.tif')).ReadAsArray()
        lgp10 = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPt10_{year}.tif')).ReadAsArray()
        lgp_equv = gdal.Open(os.path.join(work_dir, f'data_output/NB1/{country_name}_LGPEquivalent_{year}.tif')).ReadAsArray()

        
        self.applyClimaticConstraints(
            yield_input=yield_map,
            lgp=lgp,
            lgp_equv=lgp_equv,
            lgpt10=lgp10,
            omit_yld_0=True
        )
    
        # Get results
        clim_yield = self.getClimateAdjustedYield()
        fc3 = self.getClimateReductionFactor()
    
        print(f"Computed {condition_type} yield for {crop_name} using {condition_type} management")
    
        # Optional plotting
        if plot_results:
            # Load mask (assuming mask is 1 inside country, 0 outside)
            mask = gdal.Open(os.path.join(work_dir, f"data_input/{country_name}_rasterized.tif")).ReadAsArray()
            
            # Convert mask to boolean in case mask uses 255 or other codes
            mask_bool = mask.astype(bool)
            
            # Apply mask to the three layers
            yield_map_masked = np.where(mask_bool, yield_map, np.nan)
            clim_yield_masked = np.where(mask_bool, clim_yield, np.nan)
            fc3_masked = np.where(mask_bool, fc3, np.nan)
            
            # Common limits
            vmax_yield = np.nanmax([np.nanmax(clim_yield_masked), np.nanmax(yield_map_masked)])
            
            plt.figure(figsize=(22, 9))
            
            # --- Plot 1: Original Yield ---
            plt.subplot(1, 3, 1)
            plt.imshow(yield_map_masked, vmax=vmax_yield)
            plt.colorbar(shrink=0.6)
            plt.title(f"Original {condition_type.capitalize()} Yield ({crop_name})")
            
            # --- Plot 2: Climate-Constrained Yield ---
            plt.subplot(1, 3, 2)
            plt.imshow(clim_yield_masked, vmax=vmax_yield)
            plt.colorbar(shrink=0.6)
            plt.title(f"Climate-Constrained Yield ({crop_name})")
            
            # --- Plot 3: Reduction Factor Fc3 ---
            plt.subplot(1, 3, 3)
            plt.imshow(fc3_masked, vmax=1)
            plt.colorbar(shrink=0.6)
            plt.title(f"Reduction Factor Fc3 ({crop_name})")
            
            plt.tight_layout()
            plt.show()
        return clim_yield, fc3

    # Developer's Note: This code snippet below is to investigate the intermediate values used in Module III.
    #                   Do not remove this code part.

    # def getintermediate(self, i, j, yield_input, lgp, lgp_equv, lgpt10):
    #     """
    #     Generates intermediate values of Module III

    #     Returns
    #     -------
    #     TYPE : a python list.
    #         [].

    #     """
    #     lgp_agc = self.calculateLGPagc(lgp, lgp_equv)

    #     # Middle day of year for each agro-climatic constraints (used for linear interpolation purposes)
    #     mid_doy = np.array([0, 15,  45,  75, 105, 135, 165, 195, 225, 255, 285, 315, 345, 365]) # total 14 interval points


    #     if lgp_agc >=365 and self.months_P_gte_eto[i,j] == 12:
    #         gte20 = (self.gte20.drop(columns = ['365-', 'type'])).to_numpy()
    #         lt10 = (self.lt10.drop(columns=['365-', 'type'])).to_numpy()
    #         test = '365+'

    #     else:
    #         gte20 = (self.gte20.drop(columns = ['365+', 'type'])).to_numpy()
    #         lt10 = (self.lt10.drop(columns=['365+', 'type'])).to_numpy()
    #         test = '365-'
        
        
    #     # Annual mean temperature will select the relevant look-up table
        
    #     # Case I: ann_mean >= 20
    #     if self.min_T[i,j] >= 20:
    #         B_row = 1 - (np.append(0, gte20[0,:]) / 100)
    #         C_row = 1 - (np.append(0, gte20[1,:])/100)
    #         D_row = 1 - (np.append(0, gte20[2,:])/100)
        
    #     # Case II: ann_mean <= 10:
    #     elif self.min_T[i,j] <= 10:
    #         B_row = 1 - (np.append(0, lt10[0,:])/100)
    #         C_row = 1 - (np.append(0, lt10[1,:])/100)
    #         D_row = 1 - (np.append(0, lt10[2,:])/100)
        
    #     # Case III: ann_mean between 10 and 20. Linear interpolation is applied.
    #     else:

    #         # 'B' constraint row interpolation
    #         B_row_10 = np.append(0, lt10[0,:])
    #         B_row_20 = np.append(0, gte20[0,:])
    #         B_row = np.zeros(B_row_10.shape[0])

    #         for e in range(B_row_10.shape[0]):
    #             B_row[e] = 1 - ((np.interp(self.min_T[i,j], [10,20], [B_row_10[e], B_row_20[e]]))/ 100)
            
            
    #         # 'C' constraint row interpolation
    #         C_row_10 = np.append(0, lt10[1,:])
    #         C_row_20 = np.append(0, gte20[1,:])
    #         C_row = np.zeros(B_row_10.shape[0])

    #         for e in range(C_row_10.shape[0]):
    #             C_row[e] = (1 - (np.interp(self.min_T[i,j], [10,20], [C_row_10[e], C_row_20[e]]))/100)
            

    #         # 'D' constraint row interpolation
    #         D_row_10 = np.append(0, lt10[2,:])
    #         D_row_20 = np.append(0, gte20[2,:])
    #         D_row = np.zeros(B_row_10.shape[0])

    #         for e in range(D_row_10.shape[0]):
    #             D_row[e] = 1 - ((np.interp(self.min_T[i,j], [10,20], [D_row_10[e], D_row_20[e]]))/100)
        
        
    #     # 'E' constraint row interpolation
    #     E_row = np.append(0, self.lgpt10.drop(columns= 'type').iloc[0].to_numpy())
    #     E_row = 1 - (E_row/100)

    #     # Start calculation of agro-climatic constraints
    #     # 1: find agro-climatic factors of its corresponding interval of wetness days 
    #     # 2: select the most limiting factor amongst 'b', 'c', 'd' and 'e' constraints

    #     B = np.interp(lgp_agc , mid_doy, B_row)
    #     C = np.interp(lgp_agc , mid_doy, C_row)
    #     D = np.interp(lgp_agc , mid_doy, D_row)
    #     E = np.interp(lgp_agc , mid_doy, E_row)

    #     fc3 = np.round(np.min([B*C*D, E]), 2)

    #     adj_yld  = int(np.round(yield_input * fc3, 0))

    #     return [self.latitude[i,j], self.elevation[i,j], self.months_P_gte_eto[i,j], self.min_T[i,j], test, B, C, D, E, fc3, adj_yld, mid_doy, B_row, C_row, D_row, E_row, lgp_agc]

    
#----------------- End of file -------------------------#