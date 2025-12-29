"""
PyAEZ Version 3.0 (October 2025)

The `ClimateRegime` class is responsible for reading, loading, and calculating 
the agro-climatic indicators required to run PyAEZ.

Authors and Contributors:
- 2021: N. Lakmal Deshapriya
- 2022–2023: Swun Wunna Htet and Kittiphon Boonma
- 2024 (April): Swun Wunna Htet (up to version 2.3)
- 2025 (October): Dario Spiller

Modification History:

Up to Version 2.3:
1. Removed object class declarations from other modules.
2. Integrated reference water balance calculations into the routine.
3. Added new agro-climatic indicator functions.
4. Enhanced handling of monthly vs. daily time dimensions.

From Version 3.0:
1. Improved leap year management.
2. Reviewed and refined water balance calculations.
"""

import numpy as np
from pyaez.UtilitiesCalc import generateLatitudeMap, interpMonthlyToDaily, averageDailyToMonthly
from pyaez.ETOCalc import calculateETONumba, calculateNetRadiationFlux
from pyaez.LGPCalc import psh, RefWaterBalanceCalc, rainPeak, islgpt, val10day, search_cycles
from pyaez.ThermalScreening import getTempTrend, getSmoothTemp, getTemperatureGrowingPeriod
from typing import List, Tuple


np.seterr(divide='ignore', invalid='ignore') # ignore "divide by zero" or "divide by NaN" warning
np.round_ = np.round

# Initiate ClimateRegime Class instance
class ClimateRegime(object):

    def __init__(self):
        """
        Initialize the ClimateRegime object with default flag settings.
    
        Attributes:
            set_mask (bool): Flag indicating whether a spatial mask has been applied.
            set_monthly (bool): Flag indicating whether the time dimension is monthly.
            leap_year (bool): Flag indicating whether the current year is a leap year.
        """
        self.set_mask = False
        self.set_monthly = False
        self.leap_year = False
    
    def setLocationTerrainData(self, lat_min, lat_max, elevation):
        """
        (MANDATORY FUNCTION) Load geographical extents and elevation data into the class,
        and generate a latitude map based on the provided spatial dimensions.
    
        Args:
            lat_min (float): Minimum latitude of the area of interest (AOI), in decimal degrees.
            lat_max (float): Maximum latitude of the AOI, in decimal degrees.
            elevation (np.ndarray): 2D array representing elevation data in meters.
    
        Returns:
            None
        """
        self.elevation = elevation
        self.im_height = elevation.shape[0]
        self.im_width = elevation.shape[1]
        self.latitude = generateLatitudeMap(lat_min, lat_max, self.im_height, self.im_width)
        
        
    
    def setStudyAreaMask(self, admin_mask, no_data_value):
        """
        (OPTIONAL FUNCTION) Set the clipping mask for the area of interest (AOI).
    
        Args:
            admin_mask (np.ndarray): 2D binary mask used to isolate the region of interest.
            no_data_value (int): Pixel value to be treated as 'no data' and excluded from PyAEZ calculations.
    
        Returns:
            None
        """
        self.im_mask = admin_mask
        self.nodata_val = no_data_value
        self.set_mask = True

  

    def setClimateAndSoilWaterData(self, min_temp, max_temp, precipitation, short_rad, wind_speed, rel_humidity, 
                                   Sa = 100., D = 1., itflg = 1):
        """
        (MANDATORY FUNCTION) Load MONTHLY or DAILY climate data into the class and calculate:
        
        - Reference Evapotranspiration (ETo)
        - Water balance components to estimate:
            - Maximum evapotranspiration (ETm)
            - Actual evapotranspiration (ETa)
        - Temperature-based indicators for agro-climatic analysis
        
        Args:
            min_temp (3D NumPy Array): Minimum temperature [°C]
            max_temp (3D NumPy Array): Maximum temperature [°C]
            precipitation (3D NumPy Array): Total precipitation [mm/day]
            short_rad (3D NumPy Array): Solar radiation [W/m²]
            wind_speed (3D NumPy Array): Wind speed at 2m altitude [m/s]
            rel_humidity (3D NumPy Array): Relative humidity [fractional, range 0–1]
            Sa (int, float, or 2D NumPy Array, optional): Soil water holding capacity [mm/m]. Default is 100 mm.
            D (int or float, optional): Rooting depth [m]. Default is 1 m.
            itflg (int, optional): number of iterations to be performed for water balance calculation to achieve stable evaluation
        
        Returns:
            None
        """
        
        # Sanitize input data ranges
        rel_humidity[rel_humidity > 0.99] = 0.99
        rel_humidity[rel_humidity < 0.05] = 0.05
        short_rad[short_rad < 0] = 0
        wind_speed[wind_speed < 0] = 0
        
        # Time dimension check
        doy = None
        time_shapes = [min_temp, max_temp, wind_speed, short_rad, rel_humidity, precipitation]
        
        if all(arr.shape[2] == 12 for arr in time_shapes):
            self.set_monthly = True
            self.leap_year = False
            doy = 365
        elif all(arr.shape[2] == 365 for arr in time_shapes):
            self.leap_year = False
            doy = 365
        elif all(arr.shape[2] == 366 for arr in time_shapes):
            self.leap_year = True
            doy = 366
        else:
            raise ValueError("Time dimension of climate data must be 12, 365, or 366. Please check your input data.")
        
        # Initialize daily data arrays
        self.meanT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.totalPrec_daily = np.zeros((self.im_height, self.im_width, doy))
        self.minT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.maxT_daily = np.zeros((self.im_height, self.im_width, doy))
        self.shortrad_daily = np.zeros((self.im_height, self.im_width, doy))
        self.wind_daily = np.zeros((self.im_height, self.im_width, doy))
        self.rel_humidity_daily = np.zeros((self.im_height, self.im_width, doy))
        self.pet_daily = np.zeros((self.im_height, self.im_width, doy))
        self.shortrad_daily_MJm2day = np.zeros((self.im_height, self.im_width, doy))
        
        # Monthly mean temperature (used for interpolation if needed)
        meanT_monthly = (min_temp + max_temp) / 2

        if self.set_mask:
            valid_mask = self.im_mask != self.nodata_val
        else:
            valid_mask = np.ones((self.im_height, self.im_width), dtype=bool)
            
        # Use np.argwhere to iterate only over valid pixels
        for i_row, i_col in np.argwhere(valid_mask):
            if self.set_mask:
                if self.im_mask[i_row, i_col] == self.nodata_val:
                    continue
            
            if self.set_monthly:
                self.meanT_daily[i_row, i_col, :] = interpMonthlyToDaily(meanT_monthly[i_row, i_col,:], 1, doy)
                self.totalPrec_daily[i_row, i_col, :] = interpMonthlyToDaily(precipitation[i_row, i_col,:], 1, doy, no_minus_values=True)
                self.minT_daily[i_row, i_col, :] = interpMonthlyToDaily(min_temp[i_row, i_col,:], 1, doy)
                self.maxT_daily[i_row, i_col, :] = interpMonthlyToDaily(max_temp[i_row, i_col,:], 1, doy)
                self.shortrad_daily[i_row, i_col, :] = interpMonthlyToDaily(short_rad[i_row, i_col,:], 1, doy, no_minus_values=True)
                self.wind_daily[i_row, i_col, :] = interpMonthlyToDaily(wind_speed[i_row, i_col,:], 1, doy, no_minus_values=True)
                self.rel_humidity_daily[i_row, i_col, :] = interpMonthlyToDaily(rel_humidity[i_row, i_col,:], 1, doy, no_minus_values=True)
            else:
                self.meanT_daily[i_row, i_col, :] = (min_temp[i_row, i_col, :]+ max_temp[i_row, i_col, :])/2
                self.totalPrec_daily[i_row, i_col, :] = precipitation[i_row, i_col,:]
                self.minT_daily[i_row, i_col, :] = min_temp[i_row, i_col,:]
                self.maxT_daily[i_row, i_col, :] = max_temp[i_row, i_col,:]
                self.shortrad_daily[i_row, i_col, :] = short_rad[i_row, i_col,:]
                self.wind_daily[i_row, i_col, :] = wind_speed[i_row, i_col, :]
                self.rel_humidity_daily[i_row, i_col, :] = rel_humidity[i_row, i_col, :]

            # Convert radiation from W/m² to MJ/m²/day
            self.shortrad_daily_MJm2day[i_row, i_col, :] = (self.shortrad_daily[i_row, i_col, :]*3600*24)/1000000 # convert w/m2 to MJ/m2/day
            
            # Calculate reference evapotranspiration (ETo)
            self.pet_daily[i_row, i_col, :] = calculateETONumba(
                1, doy,
                self.latitude[i_row, i_col],
                self.elevation[i_row, i_col],
                self.minT_daily[i_row, i_col, :],
                self.maxT_daily[i_row, i_col, :],
                self.wind_daily[i_row, i_col, :],
                self.shortrad_daily_MJm2day[i_row, i_col, :],
                self.rel_humidity_daily[i_row, i_col, :],
                self.leap_year
            )
                
        # Sea-level adjusted mean temperature (lapse rate correction: +0.55°C per 100m elevation)
        elevation_adjustment = (self.elevation / 100) * 0.55
        self.meanT_daily_sealevel = self.meanT_daily + elevation_adjustment[:, :, np.newaxis]
        
        # Precipitation over PET ratio (avoiding division by zero and replacing NaNs with 0)
        self.P_by_PET_daily = np.divide(
            self.totalPrec_daily,
            self.pet_daily,
            out=np.zeros_like(self.pet_daily),
            where=self.pet_daily > 0
        )
        
        # Constants for snowmelt and crop coefficient - Constants for snowmelt and crop coefficient
        kc_list = np.array([0.0, 0.1, 0.2, 0.5, 1.0])  # Kc values for the reference crop
        Txsnm = 0.0   # Snow melt temperature threshold (°C)
        Fsnm = 5.5    # Snow melting coefficient
        
        # Variables initialization - Aliases (reduce attribute lookups)
        Tx365 = self.maxT_daily  # shape: (H, W, D)
        Ta365 = self.meanT_daily
        Pcp365 = self.totalPrec_daily
        self.Eto365 = self.pet_daily  # Eto
        Eto365  = self.Eto365
        
        # --- Shapes & pre-allocation ---
        H, W, T = Tx365.shape
        self.Etm365 = np.zeros((H, W, T), dtype=np.float64)
        self.Eta365 = np.zeros((H, W, T), dtype=np.float64)
        self.Sb365  = np.zeros((H, W, T), dtype=np.float64)
        self.Wb365  = np.zeros((H, W, T), dtype=np.float64)
        self.Wx365  = np.zeros((H, W, T), dtype=np.float64)
        self.kc365  = np.zeros((H, W, T), dtype=np.float64)

        # --- Valid pixel mask ---
        if self.set_mask:
            mask = (self.im_mask != self.nodata_val)
        else:
            mask = np.ones((H, W), dtype=bool)
        
        rows, cols = np.where(mask)

        #Wb_old = 0.0
        #Sb_old = 0.0
        
        # --- Main loop over valid pixels ---
        for i_row, i_col in zip(rows, cols):

            # Reset state per pixel  
            Wb_old = 0.0
            Sb_old = 0.0

            # Take 1D views to reduce indexing cost
            Tx_p  = Tx365[i_row, i_col, :].astype(np.float64, copy=False)
            Ta_p  = Ta365[i_row, i_col, :].astype(np.float64, copy=False)
            P_p   = Pcp365[i_row, i_col, :].astype(np.float64, copy=False)
            Eto_p = Eto365[i_row, i_col, :].astype(np.float64, copy=False)

            # Growing season helpers
            lgpt5_point = np.count_nonzero(Ta_p >= 5.0)

            # Determine growing season start and end based on temperature trends
            istart0, istart1 = rainPeak(Ta_p, lgpt5_point)
            
            # Get temperature trend for crop development stages
            istup = getTempTrend(Ta_p)
            
            if istup is None or len(istup) < T:
                    # Fallback / pad or recompute to match T
                    raise ValueError("getTempTrend returned invalid length for pixel ({}, {})".format(i_row, i_col))

            # --- initialize storages before spin-up
            Wb_old = 0.0
            Sb_old = 0.0

            # --- spin-up iterations
            for jj in range(itflg):  
                for t in range(T):
                    p = psh(0., self.Eto365[i_row, i_col, t])

                    Eta_new, Etm_new, Wb_new, Wx_new, Sb_new, kc_new = RefWaterBalanceCalc(
                        Tx_p[t], Ta_p[t], P_p[t],
                        Txsnm, Fsnm, Eto_p[t],
                        Wb_old, Sb_old,
                        t, istart0, istart1,
                        Sa, D, p, lgpt5_point, istup[t]
                    )

                    # Physical guards (optional)
                    Eta_new = max(Eta_new, 0.0)

                    # Update storages
                    Wb_old = Wb_new
                    Sb_old = Sb_new

                # End of 1 spin-up iteration → continue looping to stabilize
                # No outputs stored yet (unless last iteration)

            # --- After spin-up convergence (jj = itflg - 1)
            # Run one final daily loop to record results
            for t in range(T):
                p = psh(0., self.Eto365[i_row, i_col, t])

                Eta_new, Etm_new, Wb_new, Wx_new, Sb_new, kc_new = RefWaterBalanceCalc(
                    Tx_p[t], Ta_p[t], P_p[t],
                    Txsnm, Fsnm, Eto_p[t],
                    Wb_old, Sb_old,
                    t, istart0, istart1,
                    Sa, D, p, lgpt5_point, istup[t]
                )

                self.Eta365[i_row, i_col, t] = Eta_new
                self.Etm365[i_row, i_col, t] = Etm_new
                self.Wb365[i_row, i_col, t]  = Wb_new
                self.Wx365[i_row, i_col, t]  = Wx_new
                self.Sb365[i_row, i_col, t]  = Sb_new
                self.kc365[i_row, i_col, t]  = kc_new

                Wb_old = Wb_new
                Sb_old = Sb_new
    
        """
        FORTRAN CODE

        
        c calculate daily balance for water and snow bucket
        Wb365(MD365) = Wbstart
        Sb365(MD365) = Sbstart
        
        do 60 i=1,MD365
          p = psh(0, Et365(i))
          Wx365(i) = 0
        
          ! Snow period: Tmax <= Txsnm
          if (Tx365(i) .le. Txsnm) then
            kc365(i) = kc1
            etm = kc1 * Et365(i)
            Etm365(i) = etm
            sbx = sb + Pcp365(i)
            if (sbx .ge. etm) then
              Sb365(i) = sbx - etm
              Eta365(i) = etm
            else
              Sb365(i) = 0
              Eta365(i) = eta(wb, wx, etm-sbx, Sa, D, p, 0.) + sbx
            endif
            Wb365(i) = wb
            sb = Sb365(i)
        
          ! Cold water period: Ta <= 0
          else if (Ta365(i) .le. 0) then
            kc365(i) = kc2
            etm = kc2 * Et365(i)
            Etm365(i) = etm
            snm = amin1(Fsnm*(Tx365(i)-Txsnm), sb)
            sb = sb - snm
            wb = wb + snm
            sbx = sb
            if (sbx .ge. etm) then
              Sb365(i) = sbx - etm
              Eta365(i) = etm
              if (wb .gt. Sa) then
                wx = wb - Sa + Pcp365(i)
                wb = Sa
              else
                wx = Pcp365(i)
              endif
            else
              Sb365(i) = 0
              Eta365(i) = eta(wb, wx, etm-sbx, Sa, D, p, Pcp365(i)) + sbx
            endif
            Wb365(i) = wb
            sb = Sb365(i)
            Wx365(i) = wx
        
          ! Transition period: 0 < Ta < 5
          else if (Ta365(i) .lt. 5) then
            if (fromt0(i) .eq. 1) then
              kc = kc3
            else
              kc = kc7
            endif
            kc365(i) = kc
            etm = kc * Et365(i)
            Etm365(i) = etm
            snm = amin1(Fsnm*(Tx365(i)-Txsnm), sb)
            wb = wb + snm
            sb = sb - snm
            Eta365(i) = eta(wb, wx, etm, Sa, D, p, Pcp365(i))
            Sb365(i) = sb
            Wx365(i) = wx
            Wb365(i) = wb
        
          ! Warm period: Ta >= 5
          else if (Ta365(i) .ge. 5) then
            kc = kc5
            kc365(i) = kc
            etm = kc * Et365(i)
            Etm365(i) = etm
            snm = amin1(Fsnm*(Tx365(i)-Txsnm), sb)
            wb = wb + snm
            sb = sb - snm
            Eta365(i) = eta(wb, wx, etm, Sa, D, p, Pcp365(i))
            Sb365(i) = sb
            Wx365(i) = wx
            Wb365(i) = wb
          endif
        60 continue
        

        """
    def getThermalClimate(self):
        """
        Classifies rainfall and temperature seasonality into thermal climate classes.

        Returns:
            np.ndarray: 2D array of thermal climate classification codes.

        Climate Codes:
            1  - Tropical lowland (also used for equatorial uniform rainfall)
            2  - Tropical highland
            3  - Subtropical summer rainfall (also used for equatorial seasonal rainfall)
            4  - Subtropical winter rainfall
            5  - Subtropical low rainfall
            6  - Temperate oceanic
            7  - Temperate sub-continental
            8  - Temperate continental
            9  - Boreal oceanic
            10 - Boreal sub-continental
            11 - Boreal continental
            12 - Arctic

        Notes:
            - Handles both hemispheres based on latitude.
            - Special logic for equatorial zone (|latitude| ≤ 5°), mapped to existing classes.
        """
        thermal_climate = np.zeros((self.im_height, self.im_width), dtype=np.int8)

        for i_row in range(self.im_height):
            for i_col in range(self.im_width):

                # Skip masked pixels
                if self.set_mask and self.im_mask[i_row, i_col] == self.nodata_val:
                    continue

                lat = self.latitude[i_row, i_col]

                # Convert daily to monthly values
                meanT_monthly_sealevel = averageDailyToMonthly(
                    self.meanT_daily_sealevel[i_row, i_col, :], self.leap_year)
                meanT_monthly = averageDailyToMonthly(
                    self.meanT_daily[i_row, i_col, :], self.leap_year)
                P_by_PET_monthly = averageDailyToMonthly(
                    self.P_by_PET_daily[i_row, i_col, :], self.leap_year)
                monthly_precip = averageDailyToMonthly(
                    self.totalPrec_daily[i_row, i_col, :], self.leap_year)

                Ta_diff = np.max(meanT_monthly_sealevel) - np.min(meanT_monthly_sealevel)

                # Equatorial zone logic
                if abs(lat) <= 5:
                    rainfall_range = np.max(monthly_precip) - np.min(monthly_precip)
                    if rainfall_range < 50:  # Threshold for uniform rainfall
                        thermal_climate[i_row, i_col] = 1  # Tropical lowland (proxy for equatorial uniform)
                    else:
                        thermal_climate[i_row, i_col] = 3  # Subtropical summer rainfall (proxy for equatorial seasonal)
                    continue

                # Define seasonal months based on hemisphere
                if lat > 5:
                    summer_months = [3, 4, 5, 6, 7, 8]  # Apr–Sep
                    winter_months = [9, 10, 11, 0, 1, 2]  # Oct–Mar
                else:
                    summer_months = [9, 10, 11, 0, 1, 2]  # Oct–Mar
                    winter_months = [3, 4, 5, 6, 7, 8]  # Apr–Sep

                summer_PET0 = np.sum(P_by_PET_monthly[summer_months])
                winter_PET0 = np.sum(P_by_PET_monthly[winter_months])

                # Tropical climates
                if np.min(meanT_monthly_sealevel) >= 18. and Ta_diff < 15.:
                    if np.mean(meanT_monthly) < 20.:
                        thermal_climate[i_row, i_col] = 2  # Tropical highland
                    else:
                        thermal_climate[i_row, i_col] = 1  # Tropical lowland

                # Subtropical climates
                elif np.min(meanT_monthly_sealevel) >= 5. and np.sum(meanT_monthly_sealevel >= 10) >= 8:
                    total_precip = np.sum(self.totalPrec_daily[i_row, i_col, :])
                    if total_precip < 250:
                        thermal_climate[i_row, i_col] = 5  # Low rainfall
                    else:
                        summer_dominant = summer_PET0 >= winter_PET0
                        if summer_dominant:
                            thermal_climate[i_row, i_col] = 3  # Summer rainfall
                        else:
                            thermal_climate[i_row, i_col] = 4  # Winter rainfall

                # Temperate climates
                elif np.sum(meanT_monthly_sealevel >= 10) >= 4:
                    if Ta_diff <= 20:
                        thermal_climate[i_row, i_col] = 6  # Oceanic
                    elif Ta_diff <= 35:
                        thermal_climate[i_row, i_col] = 7  # Sub-continental
                    else:
                        thermal_climate[i_row, i_col] = 8  # Continental

                # Boreal climates
                elif np.sum(meanT_monthly_sealevel >= 10) >= 1:
                    if Ta_diff <= 20:
                        thermal_climate[i_row, i_col] = 9  # Oceanic
                    elif Ta_diff <= 35:
                        thermal_climate[i_row, i_col] = 10  # Sub-continental
                    else:
                        thermal_climate[i_row, i_col] = 11  # Continental

                # Arctic climate
                else:
                    thermal_climate[i_row, i_col] = 12

        return np.where(self.im_mask, thermal_climate, np.nan) if self.set_mask else thermal_climate

    
    def getThermalZone(self):
        """
        Classifies thermal zones based on monthly temperature regimes.

        Returns:
            np.ndarray: 2D array of thermal zone classification codes.

        Thermal Zone Codes:
            1  - Tropics, warm
            2  - Tropics, cool/cold/very cold
            3  - Subtropics, warm/moderately cool
            4  - Subtropics, cool
            5  - Subtropics, cold
            6  - Subtropics, very cold
            7  - Temperate, cool
            8  - Temperate, cold
            9  - Temperate, very cold
            10 - Boreal, cold
            11 - Boreal, very cold
            12 - Arctic

        Notes:
            - Classification is based on sea-level temperature thresholds and actual temperature variability.
            - Masked pixels are excluded from classification.
        """
 
        thermal_zone = np.zeros((self.im_height, self.im_width))
    
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):

                # Skip masked pixels
                if self.set_mask and self.im_mask[i_row, i_col] == self.nodata_val:
                    continue
    
                # Convert daily to monthly temperature
                meanT_monthly = averageDailyToMonthly(self.meanT_daily[i_row, i_col, :], self.leap_year)
                meanT_monthly_sealevel =  averageDailyToMonthly(self.meanT_daily_sealevel[i_row, i_col, :], self.leap_year)
    
                # Tropics
                if np.min(meanT_monthly_sealevel) >= 18 and np.max(meanT_monthly)-np.min(meanT_monthly) < 15:
                    if np.mean(meanT_monthly) > 20:
                        thermal_zone[i_row,i_col] = 1 # Tropics Warm
                    else:
                        thermal_zone[i_row,i_col] = 2 # Tropics cool/cold/very cold
                
                # Subtropics
                elif np.min(meanT_monthly_sealevel) > 5 and np.sum(meanT_monthly_sealevel > 10) >= 8:
                    if np.sum(meanT_monthly<5) >= 1 and np.sum(meanT_monthly>10) >= 4:
                        thermal_zone[i_row,i_col] =  4 # Subtropics, cool
                    elif np.sum(meanT_monthly<5) >= 1 and np.sum(meanT_monthly>10) >= 1:
                        thermal_zone[i_row,i_col] =  5 # Subtropics, cold
                    elif np.sum(meanT_monthly<10) == 12:
                        thermal_zone[i_row,i_col] =  6 # Subtropics, very cold
                    else:
                        thermal_zone[i_row,i_col] =  3 # Subtropics, warm/mod. cool
    
                # Temperate
                elif np.sum(meanT_monthly_sealevel >= 10) >= 4:
                    if np.sum(meanT_monthly<5) >= 1 and np.sum(meanT_monthly>10) >= 4:
                        thermal_zone[i_row,i_col] =  7 # Temperate, cool
                    elif np.sum(meanT_monthly<5) >= 1 and np.sum(meanT_monthly>10) >= 1:
                        thermal_zone[i_row,i_col] =  8 # Temperate, cold
                    elif np.sum(meanT_monthly<10) == 12:
                        thermal_zone[i_row,i_col] =  9 # Temperate, very cold
    
                # Boreal
                elif np.sum(meanT_monthly_sealevel >= 10) >= 1:
                    if np.sum(meanT_monthly<5) >= 1 and np.sum(meanT_monthly>10) >= 1:
                        thermal_zone[i_row,i_col] = 10 # Boreal, cold
                    elif np.sum(meanT_monthly<10) == 12:
                        thermal_zone[i_row,i_col] = 11 # Boreal, very cold
                
                # Arctic
                else:
                        thermal_zone[i_row,i_col] = 12 # Arctic
    
        if self.set_mask:
            return np.where(self.im_mask, thermal_zone, np.nan)
        else:
            return thermal_zone

    def _computeThermalLGP(self, threshold, attr_name):
        """
        Internal method to compute the Thermal Length of Growing Period (LGP)
        for a given temperature threshold.

        Args:
            threshold (float): Temperature threshold in °C.
            attr_name (str): Name of the attribute to store the result (e.g., 'lgpt0').

        Returns:
            np.ndarray: 2D array of number of days with mean temperature ≥ threshold.
        """
        lgp = np.sum(self.meanT_daily >= threshold, axis=2)

        if self.set_mask:
            lgp = np.where(self.im_mask, lgp, np.nan)

        setattr(self, attr_name, lgp.copy())
        return lgp

    def getThermalLGP0(self):
        """
        Calculates the Thermal Length of Growing Period (LGP) using a temperature threshold of 0°C.
        Returns:
            np.ndarray: Days with mean temperature ≥ 0°C.
        """
        return self._computeThermalLGP(threshold=0, attr_name='lgpt0')


    def getThermalLGP5(self):
        """
        Calculates the Thermal Length of Growing Period (LGP) using a temperature threshold of 5°C.
        Returns:
            np.ndarray: Days with mean temperature ≥ 5°C.
        """
        return self._computeThermalLGP(threshold=5, attr_name='lgpt5')


    def getThermalLGP10(self):
        """
        Calculates the Thermal Length of Growing Period (LGP) using a temperature threshold of 10°C.
        Returns:
            np.ndarray: Days with mean temperature ≥ 10°C.
        """
        return self._computeThermalLGP(threshold=10, attr_name='lgpt10')

    def _computeTemperatureSum(self, threshold, attr_name):
        """
        Internal method to compute temperature summation above a given threshold.

        Args:
            threshold (float): Temperature threshold in °C.
            attr_name (str): Name of the attribute to store the result (e.g., 'tsum0').

        Returns:
            np.ndarray: 2D array of accumulated daily mean temperatures above the threshold.
                        Units: Degree-Days
        """
        tempT = self.meanT_daily.copy()
        tempT[tempT < threshold] = 0
        tsum = np.round(np.sum(tempT, axis=2), decimals=0)

        if self.set_mask:
            tsum = np.where(self.im_mask, tsum, np.nan)

        setattr(self, attr_name, tsum.copy())
        return tsum

    def getTemperatureSum0(self):
        """
        Calculates temperature summation for days with mean temperature ≥ 0°C.

        Returns:
            np.ndarray: Accumulated degree-days above 0°C.
        """
        return self._computeTemperatureSum(threshold=0, attr_name='tsum0')


    def getTemperatureSum5(self):
        """
        Calculates temperature summation for days with mean temperature ≥ 5°C.

        Returns:
            np.ndarray: Accumulated degree-days above 5°C.
        """
        return self._computeTemperatureSum(threshold=5, attr_name='tsum5')


    def getTemperatureSum10(self):
        """
        Calculates temperature summation for days with mean temperature ≥ 10°C.

        Returns:
            np.ndarray: Accumulated degree-days above 10°C.
        """
        return self._computeTemperatureSum(threshold=10, attr_name='tsum10')

    def getTemperatureProfile(self):
        """
        Classifies temperature profile based on daily temperature transitions.

        Returns:
            list of np.ndarray: 18 2D arrays representing the number of days per pixel
                                where temperature is rising (A1–A9) or falling (B1–B9)
                                within specific temperature ranges.

        Temperature Profile Classes:
            A1–A9: Warming transitions
            B1–B9: Cooling transitions

            Ranges:
                - A1/B1: ≥ 30°C
                - A2/B2: 25–30°C
                - A3/B3: 20–25°C
                - A4/B4: 15–20°C
                - A5/B5: 10–15°C
                - A6/B6: 5–10°C
                - A7/B7: 0–5°C
                - A8/B8: -5–0°C
                - A9/B9: < -5°C

        Notes:
            - Uses a 5th-degree polynomial fit to smooth daily temperature time series.
            - Handles leap years.
            - Applies masking if enabled.
        """
        def compute_transition_counts(meanT_first, meanT_diff, low, high, direction):
            """Helper to compute transition counts for warming or cooling."""
            if direction == "warming":
                mask = np.logical_and(meanT_diff > 0, np.logical_and(meanT_first >= low, meanT_first < high))
            elif direction == "cooling":
                mask = np.logical_and(meanT_diff < 0, np.logical_and(meanT_first >= low, meanT_first < high))
            else:
                raise ValueError("Direction must be 'warming' or 'cooling'")
            result = np.sum(mask, axis=2)
            if self.set_mask:
                result = np.ma.masked_where(self.im_mask == 0, result)
            return result

        DAYS_IN_YEAR = 366 if self.leap_year else 365
        days = np.arange(DAYS_IN_YEAR)

        # Interpolate daily temperature using polynomial fit
        interp_daily_temp = np.zeros((self.im_height, self.im_width, DAYS_IN_YEAR))
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):
                temp_1D = self.meanT_daily[i_row, i_col, :]
                poly_fit = np.poly1d(np.polyfit(days, temp_1D, 5))
                interp_daily_temp[i_row, i_col, :] = poly_fit(days)

        # Extend time series by one day to compute daily differences
        extended_temp = np.concatenate((interp_daily_temp, interp_daily_temp[:, :, :1]), axis=-1)
        meanT_first = extended_temp[:, :, :-1]
        meanT_diff = extended_temp[:, :, 1:] - meanT_first

        # Define temperature bins from warmest to coldest
        bins = [(30, np.inf), (25, 30), (20, 25), (15, 20), (10, 15),
                (5, 10), (0, 5), (-5, 0), (-np.inf, -5)]

        # Explicit variable names for warming transitions
        A1 = compute_transition_counts(meanT_first, meanT_diff, *bins[0], "warming")
        A2 = compute_transition_counts(meanT_first, meanT_diff, *bins[1], "warming")
        A3 = compute_transition_counts(meanT_first, meanT_diff, *bins[2], "warming")
        A4 = compute_transition_counts(meanT_first, meanT_diff, *bins[3], "warming")
        A5 = compute_transition_counts(meanT_first, meanT_diff, *bins[4], "warming")
        A6 = compute_transition_counts(meanT_first, meanT_diff, *bins[5], "warming")
        A7 = compute_transition_counts(meanT_first, meanT_diff, *bins[6], "warming")
        A8 = compute_transition_counts(meanT_first, meanT_diff, *bins[7], "warming")
        A9 = compute_transition_counts(meanT_first, meanT_diff, *bins[8], "warming")

        # Explicit variable names for cooling transitions
        B1 = compute_transition_counts(meanT_first, meanT_diff, *bins[0], "cooling")
        B2 = compute_transition_counts(meanT_first, meanT_diff, *bins[1], "cooling")
        B3 = compute_transition_counts(meanT_first, meanT_diff, *bins[2], "cooling")
        B4 = compute_transition_counts(meanT_first, meanT_diff, *bins[3], "cooling")
        B5 = compute_transition_counts(meanT_first, meanT_diff, *bins[4], "cooling")
        B6 = compute_transition_counts(meanT_first, meanT_diff, *bins[5], "cooling")
        B7 = compute_transition_counts(meanT_first, meanT_diff, *bins[6], "cooling")
        B8 = compute_transition_counts(meanT_first, meanT_diff, *bins[7], "cooling")
        B9 = compute_transition_counts(meanT_first, meanT_diff, *bins[8], "cooling")

        return [A1, A2, A3, A4, A5, A6, A7, A8, A9,
                B1, B2, B3, B4, B5, B6, B7, B8, B9]
    
    def val10day_2years(arr):
        """
        Python equivalent of the Fortran subroutine val10day.
        Computes a 10-day trailing average with wrap-around logic.
        Returns an array of twice the input length: first half with averages,
        second half as a duplicate of the first.
        """
        n = len(arr)
        arr_padded = np.concatenate((arr[:9], arr))  # pad first 9 days from start
        val10 = np.zeros(n * 2)

        # Compute 10-day trailing average
        for jd in range(9, n + 9):
            window = arr_padded[jd - 9: jd + 1]  # 10-day window ending at jd
            val10[jd] = np.mean(window)

        # Wrap-around for first 9 days
        for jd in range(9):
            val10[jd] = val10[n + jd]

        # Copy first year to second year
        val10[n:] = val10[:n]

        return val10

    def getLGP(self):
        """Calculate length of growing period (LGP).

        Args:
            None.
        Return:
           lgp (2D-NumPy Array): length of growing periods [Unit: Days].
        """        
        lgp_tot = np.zeros((self.im_height, self.im_width))
        #============================
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_row, i_col] == self.nodata_val:
                        continue
                Etm365X = np.append(self.Etm365[i_row, i_col, :], self.Etm365[i_row, i_col, :])
                Eta365X = np.append(self.Eta365[i_row, i_col, :], self.Eta365[i_row, i_col, :])
                islgp = islgpt(self.meanT_daily[i_row, i_col, :])
                xx = val10day(Eta365X)
                yy = val10day(Etm365X)
                lgp_whole = xx[:365]/yy[:365]
                count = 0
                for i in range(len(lgp_whole)):
                    if islgp[i] == 1 and lgp_whole[i] >= 0.4:
                        count = count+1

                lgp_tot[i_row, i_col] = count

        if self.set_mask:
            return np.where(self.im_mask, lgp_tot, np.nan)
        else:
            return lgp_tot
  
    
    def getLGPClassified(self, lgp): # Original PyAEZ source code
        """This function calculates the classification of moisture regimes based on LGP.

        Args:
            lgp (2D-NumPy Array): Length of Growing Period [Unit: Days]

        Return:
            lgp_class (2D-NumPy Array): Moisture regime classes.

        """        
        lgp_class = np.zeros(lgp.shape)

        lgp_class[lgp>=365] = 7 # Per-humid
        lgp_class[np.logical_and(lgp>=270, lgp<365)] = 6 # Humid
        lgp_class[np.logical_and(lgp>=180, lgp<270)] = 5 # Sub-humid
        lgp_class[np.logical_and(lgp>=120, lgp<180)] = 4 # Moist semi-arid
        lgp_class[np.logical_and(lgp>=60, lgp<120)] = 3 # Dry semi-arid
        lgp_class[np.logical_and(lgp>0, lgp<60)] = 2 # Arid
        lgp_class[lgp<=0] = 1 # Hyper-arid

        if self.set_mask:
            return np.where(self.im_mask, lgp_class, np.nan)
        else:
            return lgp_class
        
        
    def getLGPEquivalent(self): 
        """Calculate the equivalent length of growing period.

        Args:
            None.
        Return:
            lgp_equv (2D-NumPy Array): equivalent length of growing period [Unit: Days].
        """        
        moisture_index = np.sum(self.totalPrec_daily, axis=2)/np.sum(self.pet_daily, axis=2)

        lgp_equv = 14.0 + 293.66*moisture_index - 61.25*moisture_index*moisture_index
        lgp_equv[moisture_index > 2.4] = 366

        if self.set_mask:
            return np.where(self.im_mask, lgp_equv, np.nan)
        else:
            return lgp_equv

        # '''
        # Existing Issue: The moisture index calculation is technical aligned with FORTRAN routine, 
        # results are still different from GAEZ; causing large discrepancy. 
        # Overall, there are no changes with the calculation steps and logics.
        # '''
      
    def TZoneFallowRequirement(self, tzone):
        """
        The function calculates the temperature zones applied for fallow requirements which 
        requires thermal zone to classify. 

        Args:
            tzone (2D-NumPy Array): thermal zone classes.
        Return:
            tzone_fallow (2D-NumPy Array): thermal zone for fallow requirements.

        """

        # the algorithm needs to calculate the annual mean temperature.
        tzonefallow = np.zeros((self.im_height, self.im_width), dtype= int)
        annual_Tmean = np.mean(self.meanT_daily, axis = 2)

        # thermal zone class definitions for fallow requirement
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):

                if self.set_mask:
                    if self.im_mask[i_row, i_col] == self.nodata_val:
                        continue
                # Checking tropics thermal zone
                if tzone[i_row, i_col] == 1 or tzone[i_row, i_col] == 2:
                    
                    # Class 1: tropics, mean annual T > 25 deg C
                    if annual_Tmean[i_row, i_col] > 25:
                        tzonefallow[i_row, i_col] = 1
                    
                    # Class 2: tropics, mean annual T 20-25 deg C
                    elif annual_Tmean[i_row, i_col] > 20:
                        tzonefallow[i_row, i_col] = 2
                    
                    # Class 3: tropics, mean annual T 15-20 deg C
                    elif annual_Tmean[i_row, i_col] > 15:
                        tzonefallow[i_row, i_col] = 3
                    
                    # Class 4: tropics, mean annual T < 15 deg C
                    else:
                        tzonefallow[i_row, i_col] = 4
                
                # Checking the non-tropical zones
                else:
                    meanT_monthly = averageDailyToMonthly(self.meanT_daily[i_row, i_col, :], self.leap_year)
                    # Class 5: mean T of the warmest month > 20 deg C
                    if np.max(meanT_monthly) > 20:
                        tzonefallow[i_row, i_col] = 5
                        
                    else:
                        tzonefallow[i_row, i_col] = 6
                            
        if self.set_mask:
            return np.where(self.im_mask, tzonefallow, np.nan)
        else:
            return tzonefallow
    
    def AirFrostIndexandPermafrostEvaluation(self):
        """
        The function calculates the air frost index which is used for evaluation of 
        occurrence of continuous or discontinuous permafrost condtions executed in 
        GAEZ v4. Two outputs of numerical air frost index and classified reference
        permafrost zones are returned.

        Args:
            None.
        Return:
            air_frost_index/permafrost : a python list: [air frost number, permafrost classes]

        """
        fi = np.zeros((self.im_height, self.im_width), dtype=float)
        permafrost = np.zeros((self.im_height, self.im_width), dtype=int)
        ddt = np.zeros((self.im_height, self.im_width), dtype=float) # thawing index
        ddf = np.zeros((self.im_height, self.im_width), dtype=float) # freezing index
        meanT_gt_0 = self.meanT_daily.copy()
        meanT_le_0 = self.meanT_daily.copy()
        
        meanT_gt_0[meanT_gt_0 <=0] = 0 # removing all negative temperatures for summation
        meanT_le_0[meanT_gt_0 >0] = 0 # removing all positive temperatures for summation 
        ddt = np.sum(meanT_gt_0, axis = 2)
        ddf = - np.sum(meanT_le_0, axis = 2)  
        fi = np.sqrt(ddf)/(np.sqrt(ddf) + np.sqrt(ddt)) 
        # now, we will classify the permafrost zones (Reference: GAEZ v4 model documentation: Pg35 -37)
        for i_row in range(self.im_height):
            for i_col in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_row, i_col] == self.nodata_val:
                        continue         
                # Continuous Permafrost Class
                if fi[i_row, i_col]> 0.625:
                    permafrost[i_row, i_col] = 1
                
                # Discontinuous Permafrost Class
                if fi[i_row, i_col]> 0.57 and fi[i_row, i_col]< 0.625:
                    permafrost[i_row, i_col] = 2
                
                # Sporadic Permafrost Class
                if fi[i_row, i_col]> 0.495 and fi[i_row, i_col]< 0.57:
                    permafrost[i_row, i_col] = 3
                
                # No Permafrost Class
                if fi[i_row, i_col]< 0.495:
                    permafrost[i_row, i_col] = 4
        # to remove the division by zero, the nan values will be converted into
        fi = np.nan_to_num(fi)

        if self.set_mask:
            return [np.where(self.im_mask, fi, np.nan), np.where(self.im_mask, permafrost , np.nan)]
        else:
            return [fi, permafrost]
        
    
    def AEZClassification(self, tclimate, lgp, lgp_equv, lgpt_5, soil_terrain_lulc, permafrost):
        """The AEZ inventory combines spatial layers of thermal and moisture regimes 
        with broad categories of soil/terrain qualities.

        Args:
            tclimate (2D-NumPy Array): Thermal Climate classes
            lgp (2D-NumPy Array): Length of Growing Period [Unit: Days]
            lgp_equv (2D-NumPy Array): Equivalent length of growing periods [Unit:Days]
            lgpt_5 (2D-NumPy Array): Thermal LGP of days with Ta>5˚C
            soil_terrain_lulc (2D-NumPy Array): soil/terrain/special land cover classes (8 classes)
            permafrost (2D-NumPy Array): Permafrost classes

        Return:
            aez (2D-NumPy Array): aez zones (57 classes).
        """        
        
        #1st step: reclassifying the existing 12 classes of thermal climate into 6 major thermal climate.
        # Class 1: Tropics, lowland
        # Class 2: Tropics, highland
        # Class 3: Subtropics
        # Class 4: Temperate Climate
        # Class 5: Boreal Climate
        # Class 6: Arctic Climate
    
        aez_tclimate = np.zeros((self.im_height, self.im_width), dtype=int)

        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_r, i_c] == self.nodata_val:
                        continue

                    else:

                        # tropics highland
                        if tclimate[i_r, i_c] == 1:
                            aez_tclimate[i_r, i_c] = 1

                        elif tclimate[i_r, i_c] == 2:
                            aez_tclimate[i_r, i_c] = 2

                        elif tclimate[i_r, i_c] == 3:
                            aez_tclimate[i_r, i_c] = 3

                        elif tclimate[i_r, i_c] == 4:
                            aez_tclimate[i_r, i_c] = 3

                        elif tclimate[i_r, i_c] == 5:
                            aez_tclimate[i_r, i_c] = 3

                        # grouping all the temperate classes into a single class 4
                        elif tclimate[i_r, i_c] == 6:
                            aez_tclimate[i_r, i_c] = 4

                        elif tclimate[i_r, i_c] == 7:
                            aez_tclimate[i_r, i_c] = 4

                        elif tclimate[i_r, i_c] == 8:
                            aez_tclimate[i_r, i_c] = 4

                        # grouping all the boreal classes into a single class 5
                        elif tclimate[i_r, i_c] == 9:
                            aez_tclimate[i_r, i_c] = 5

                        elif tclimate[i_r, i_c] == 10:
                            aez_tclimate[i_r, i_c] = 5

                        elif tclimate[i_r, i_c] == 11:
                            aez_tclimate[i_r, i_c] = 5

                        # changing the arctic class into class 6
                        elif tclimate[i_r, i_c] == 12:
                            aez_tclimate[i_r, i_c] = 6

        # 2nd Step: Classification of Thermal Zones
        aez_tzone = np.zeros((self.im_height, self.im_width), dtype=int)


        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_r, i_c] == self.nodata_val:
                        continue
                    else:
                        mean_temp = np.copy(self.meanT_daily[i_r, i_c, :])
                        meanT_monthly = averageDailyToMonthly(
                            mean_temp, self.leap_year)
                        # one conditional parameter for temperature accumulation
                        temp_acc_10deg = np.copy(self.meanT_daily[i_r, i_c, :])
                        temp_acc_10deg[temp_acc_10deg < 10] = 0

                        # Warm Tzone (TZ1)
                        if np.sum(meanT_monthly >= 10) == 12 and np.mean(mean_temp) >= 20:
                            aez_tzone[i_r, i_c] = 1

                        # Moderately cool Tzone (TZ2)
                        elif np.sum(meanT_monthly >= 5) == 12 and np.sum(meanT_monthly >= 10) >= 8:
                            aez_tzone[i_r, i_c] = 2

                        # TZ3 Moderate
                        elif aez_tclimate[i_r, i_c] == 4 and np.sum(meanT_monthly >= 10) >= 5 and np.sum(mean_temp > 20) >= 75 and np.sum(temp_acc_10deg) > 3000:
                            aez_tzone[i_r, i_c] = 3

                        # TZ4 Cool
                        elif np.sum(meanT_monthly >= 10) >= 4 and np.mean(mean_temp) >= 0:
                            aez_tzone[i_r, i_c] = 4

                        # TZ5 Cold
                        elif np.sum(meanT_monthly >= 10) in range(1, 4) and np.mean(mean_temp) >= 0:
                            aez_tzone[i_r, i_c] = 5

                        # TZ6 Very cold
                        elif np.sum(meanT_monthly < 10) == 12 or np.mean(mean_temp) < 0:
                            aez_tzone[i_r, i_c] = 6

        # 3rd Step: Creation of Temperature Regime Classes
        # Temperature Regime Class Definition
        # 1 = Tropics, lowland (TRC1)
        # 2 = Tropics, highland (TRC2)
        # 3 = Subtropics, warm (TRC3)
        # 4 = Subtropics, moderately cool (TRC4)
        # 5 = Subtropics, cool (TRC5)
        # 6 = Temperate, moderate (TRC6)
        # 7 = Temperate, cool (TRC7)
        # 8 = Boreal, cold, no continuous or discontinuous occurrence of permafrost (TRC8)
        # 9 = Boreal, cold, with continuous or discontinuous occurrence of permafrost (TRC9)
        # 10 = Arctic, very cold (TRC10)

        aez_temp_regime = np.zeros((self.im_height, self.im_width), dtype=int)

        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_r, i_c] == self.nodata_val:
                        continue
                    else:

                        if aez_tclimate[i_r, i_c] == 1 and aez_tzone[i_r, i_c] == 1:
                            aez_temp_regime[i_r, i_c] = 1  # Tropics, lowland

                        elif aez_tclimate[i_r, i_c] == 2 and aez_tzone[i_r, i_c] in [2, 4]:
                            aez_temp_regime[i_r, i_c] = 2  # Tropics, highland

                        elif aez_tclimate[i_r, i_c] == 3 and aez_tzone[i_r, i_c] == 1:
                            aez_temp_regime[i_r, i_c] = 3  # Subtropics, warm

                        elif aez_tclimate[i_r, i_c] == 3 and aez_tzone[i_r, i_c] == 2:
                            # Subtropics,moderate cool
                            aez_temp_regime[i_r, i_c] = 4

                        elif aez_tclimate[i_r, i_c] == 3 and aez_tzone[i_r, i_c] == 4:
                            aez_temp_regime[i_r, i_c] = 5  # Subtropics,cool

                        elif aez_tclimate[i_r, i_c] == 4 and aez_tzone[i_r, i_c] == 3:
                            # Temperate, moderate
                            aez_temp_regime[i_r, i_c] = 6

                        elif aez_tclimate[i_r, i_c] == 4 and aez_tzone[i_r, i_c] == 4:
                            aez_temp_regime[i_r, i_c] = 7  # Temperate, cool

                        elif aez_tclimate[i_r, i_c] in range(2, 6) and aez_tzone[i_r, i_c] == 5:
                            if np.logical_or(permafrost[i_r, i_c] == 1, permafrost[i_r, i_c] == 2) == False:
                                # Boreal/Cold, no
                                aez_temp_regime[i_r, i_c] = 8
                            else:
                                # Boreal/Cold, with permafrost
                                aez_temp_regime[i_r, i_c] = 9

                        elif aez_tclimate[i_r, i_c] in range(2, 7) and aez_tzone[i_r, i_c] == 6:
                            aez_temp_regime[i_r, i_c] = 10  # Arctic/Very Cold

        # 4th Step: Moisture Regime classes
        # Moisture Regime Class Definition
        # 1 = M1 (desert/arid areas, 0 <= LGP* < 60)
        # 2 = M2 (semi-arid/dry areas, 60 <= LGP* < 180)
        # 3 = M3 (sub-humid/moist areas, 180 <= LGP* < 270)
        # 4 = M4 (humid/wet areas, LGP* >= 270)

        aez_moisture_regime = np.zeros(
            (self.im_height, self.im_width), dtype=int)

        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_r, i_c] == self.nodata_val:
                        continue
                    else:

                        # check if LGP t>5 is greater or less than 330 days. If greater, LGP will be used; otherwise, LGP_equv will be used.
                        if lgpt_5[i_r, i_c] > 330:

                            # Class 4 (M4)
                            if lgp[i_r, i_c] >= 270:
                                aez_moisture_regime[i_r, i_c] = 4

                            # Class 3 (M3)
                            elif lgp[i_r, i_c] >= 180 and lgp[i_r, i_c] < 270:
                                aez_moisture_regime[i_r, i_c] = 3

                            # Class 2 (M2)
                            elif lgp[i_r, i_c] >= 60 and lgp[i_r, i_c] < 180:
                                aez_moisture_regime[i_r, i_c] = 2

                            # Class 1 (M1)
                            elif lgp[i_r, i_c] >= 0 and lgp[i_r, i_c] < 60:
                                aez_moisture_regime[i_r, i_c] = 1

                        elif lgpt_5[i_r, i_c] <= 330:

                            # Class 4 (M4)
                            if lgp_equv[i_r, i_c] >= 270:
                                aez_moisture_regime[i_r, i_c] = 4

                            # Class 3 (M3)
                            elif lgp_equv[i_r, i_c] >= 180 and lgp_equv[i_r, i_c] < 270:
                                aez_moisture_regime[i_r, i_c] = 3

                            # Class 2 (M2)
                            elif lgp_equv[i_r, i_c] >= 60 and lgp_equv[i_r, i_c] < 180:
                                aez_moisture_regime[i_r, i_c] = 2

                            # Class 1 (M1)
                            elif lgp_equv[i_r, i_c] >= 0 and lgp_equv[i_r, i_c] < 60:
                                aez_moisture_regime[i_r, i_c] = 1

        # Now, we will classify the agro-ecological zonation
        # By GAEZ v4 Documentation, there are prioritized sequential assignment of AEZ classes in order to ensure the consistency of classification
        aez = np.zeros((self.im_height, self.im_width), dtype=int)

        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                if self.set_mask:
                    if self.im_mask[i_r, i_c] == self.nodata_val:
                        continue
                    else:
                        # if it's urban built-up lulc, Dominantly urban/built-up land
                        if soil_terrain_lulc[i_r, i_c] == 8:
                            aez[i_r, i_c] = 56

                        # if it's water/ dominantly water
                        elif soil_terrain_lulc[i_r, i_c] == 7:
                            aez[i_r, i_c] = 57

                        # if it's dominantly very steep terrain/Dominantly very steep terrain
                        elif soil_terrain_lulc[i_r, i_c] == 1:
                            aez[i_r, i_c] = 49

                        # if it's irrigated soils/ Land with ample irrigated soils
                        elif soil_terrain_lulc[i_r, i_c] == 6:
                            aez[i_r, i_c] = 51

                        # if it's hydromorphic soils/ Dominantly hydromorphic soils
                        elif soil_terrain_lulc[i_r, i_c] == 2:
                            aez[i_r, i_c] = 52

                        # Desert/Arid climate
                        elif aez_moisture_regime[i_r, i_c] == 1:
                            aez[i_r, i_c] = 53

                        # BO/Cold climate, with Permafrost
                        elif aez_temp_regime[i_r, i_c] == 9 and aez_moisture_regime[i_r, i_c] in [1, 2, 3, 4] == True:
                            aez[i_r, i_c] = 54

                        # Arctic/ Very cold climate
                        elif aez_temp_regime[i_r, i_c] == 10 and aez_moisture_regime[i_r, i_c] in [1, 2, 3, 4] == True:
                            aez[i_r, i_c] = 55

                        # Severe soil/terrain limitations
                        elif soil_terrain_lulc[i_r, i_c] == 5:
                            aez[i_r, i_c] = 50

                        #######
                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 1

                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 2

                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 3

                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 4

                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 5

                        elif aez_temp_regime[i_r, i_c] == 1 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 6
                        ####
                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 7

                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 8

                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 9

                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 10

                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 11

                        elif aez_temp_regime[i_r, i_c] == 2 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 12
                        ###
                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 13

                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 14

                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 15

                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 16

                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 17

                        elif aez_temp_regime[i_r, i_c] == 3 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 18
                        #####
                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 19

                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 20

                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 21

                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 22

                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 23

                        elif aez_temp_regime[i_r, i_c] == 4 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 24
                        #####
                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 25

                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 26

                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 27

                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 28

                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 29

                        elif aez_temp_regime[i_r, i_c] == 5 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 30
                        ######

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 31

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 32

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 33

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 34

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 35

                        elif aez_temp_regime[i_r, i_c] == 6 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 36

                        ###
                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 37

                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 38

                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 39

                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 40

                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 41

                        elif aez_temp_regime[i_r, i_c] == 7 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 42
                        #####

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 43

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 2 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 44

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 45

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 3 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 46

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 3:
                            aez[i_r, i_c] = 47

                        elif aez_temp_regime[i_r, i_c] == 8 and aez_moisture_regime[i_r, i_c] == 4 and soil_terrain_lulc[i_r, i_c] == 4:
                            aez[i_r, i_c] = 48          

        if self.set_mask:
            return np.where(self.im_mask, aez, np.nan)
        else:        
            return aez
    
    """ 
    Note from Swun: In this code, the logic of temperature amplitude is not added 
    as it brings big discrepency in the temperature regime calculation (in India) 
    compared to previous code. However, the classification schema is now adjusted 
    according to Gunther's agreement and the documentation.
    """
         
    def getMultiCroppingZones(self, t_climate, lgp, lgp_t5, lgp_t10, ts_t10, ts_t0):
        
        # defining the constant arrays for rainfed and irrigated conditions, all pixel values start with 1
        multi_crop_rain = np.zeros((self.im_height, self.im_width), dtype = int) # all values started with Zone A
        multi_crop_irr = np.zeros((self.im_height, self.im_width), dtype = int) # all vauels starts with Zone A
        
        ts_g_t5 = np.zeros((self.im_height, self.im_width))
        ts_g_t10 = np.zeros((self.im_height, self.im_width))

        if self.leap_year:
            DAYS_IN_YEAR = 366
        else:
            DAYS_IN_YEAR = 365
            
        # Calculation of Accumulated temperature during the growing period at specific temperature thresholds: 5 and 10 degree Celsius
        
        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                
                if self.set_mask:
                    
                    if self.im_mask[i_r, i_c]== self.nodata_val:
                        continue
                    
                    else:
                        
                        temp_1D = self.meanT_daily[i_r, i_c, :]
                        days = np.arange(0,DAYS_IN_YEAR)
                        
                        deg = 5 # order of polynomical fit
                        
                        # creating the function of polyfit
                        polyfit = np.poly1d(np.polyfit(days,temp_1D,deg))
                        
                        # getting the interpolated value at each DOY
                        interp_daily_temp = polyfit(days)
                        
                        # Getting the start and end day of vegetative period
                        # The crop growth requires minimum temperature of at least 5 deg Celsius
                        # If not, the first DOY and the lst DOY of a year will be considered
                        try:
                            veg_period = days[interp_daily_temp >=5]
                            start_veg = veg_period[0]
                            end_veg = veg_period[-1]
                        except:
                            start_veg = 0
                            end_veg = DAYS_IN_YEAR - 1
                        
                        # Slicing the temperature within the vegetative period
                        interp_meanT_veg_T5 = interp_daily_temp[start_veg:end_veg]
                        interp_meanT_veg_T10 =  interp_daily_temp[start_veg:end_veg] *1
                        
                        # Removing the temperature of 5 and 10 deg Celsius thresholds
                        interp_meanT_veg_T5[interp_meanT_veg_T5 < 5] = 0
                        interp_meanT_veg_T10[interp_meanT_veg_T10 <10] = 0
                        
                        # Calculation of Accumulated temperatures during growing period
                        ts_g_t5[i_r, i_c] = np.sum(interp_meanT_veg_T5)
                        ts_g_t10[i_r, i_c] = np.sum(interp_meanT_veg_T10)
        
        """Multi cropping zonation for rainfed conditions"""
        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                
                if self.set_mask:
                    
                    if self.im_mask[i_r, i_c]== self.nodata_val:
                        continue
                    
                    else:
                        
                        if t_climate[i_r, i_c]== 1:
                            
                            if np.all([lgp[i_r, i_c]>=360, lgp_t5[i_r, i_c]>=360, lgp_t10[i_r, i_c]>=360, ts_t0[i_r, i_c]>=7200, ts_t10[i_r, i_c]>=7000])== True:
                                multi_crop_rain[i_r, i_c] = 8
                            
                            elif np.all([lgp[i_r, i_c]>=300, lgp_t5[i_r, i_c]>=300, lgp_t10[i_r, i_c]>=240, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=5100, ts_g_t10[i_r, i_c]>=4800])== True:
                                multi_crop_rain[i_r, i_c] = 6
                            
                            elif np.all([lgp[i_r, i_c]>=270, lgp_t5[i_r, i_c]>=270, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=5500, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_rain[i_r, i_c] = 4 # Ok
                                
                            elif np.all([lgp[i_r, i_c]>=240, lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=6400, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_rain[i_r, i_c] = 4 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=210, lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_rain[i_r, i_c] = 4 # OK
                            
                            elif np.all([lgp[i_r, i_c]>=220, lgp_t5[i_r, i_c]>=220, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=5500, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])== True:
                                multi_crop_rain[i_r, i_c] = 3 #OK
                            
                            elif np.all([lgp[i_r, i_c]>=200, lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=6400, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])== True:
                                multi_crop_rain[i_r, i_c] = 3# OK
                            
                            elif np.all([lgp[i_r, i_c]>=180, lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])== True:
                                multi_crop_rain[i_r, i_c] = 3 # OK
                            
                            elif np.all([lgp[i_r, i_c]>=45, lgp_t5[i_r, i_c]>=120, lgp_t10[i_r, i_c]>=90, ts_t0[i_r, i_c]>=1600, ts_t10[i_r, i_c]>=1200]) == True:
                                multi_crop_rain[i_r, i_c] = 2 # Ok
                                
                            else:
                                multi_crop_rain[i_r, i_c] = 1 # Ok
                            
                        elif t_climate[i_r, i_c] != 1:
                            
                            if np.all([lgp[i_r, i_c]>=360, lgp_t5[i_r, i_c]>=360, lgp_t10[i_r, i_c]>=330, ts_t0[i_r, i_c]>=7200, ts_t10[i_r, i_c]>=7000])== True:
                                multi_crop_rain[i_r, i_c] = 8 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=330, lgp_t5[i_r, i_c]>=330, lgp_t10[i_r, i_c]>=270, ts_t0[i_r, i_c]>=5700, ts_t10[i_r, i_c]>=5500])== True:
                                multi_crop_rain[i_r, i_c] = 7 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=300, lgp_t5[i_r, i_c]>=300, lgp_t10[i_r, i_c]>=240, ts_t0[i_r, i_c]>=5400, ts_t10[i_r, i_c]>=5100, ts_g_t5[i_r, i_c]>=5100, ts_g_t10[i_r, i_c]>=4800])== True:
                                multi_crop_rain[i_r, i_c] = 6 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=240, lgp_t5[i_r, i_c]>=270, lgp_t10[i_r, i_c]>=180, ts_t0[i_r, i_c]>=4800, ts_t10[i_r, i_c]>=4500, ts_g_t5[i_r, i_c]>=4300, ts_g_t10[i_r, i_c]>=4000])== True:
                                multi_crop_rain[i_r, i_c] = 5 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=210, lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=4500, ts_t10[i_r, i_c]>=3600, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_rain[i_r, i_c] = 4 #OK
                            
                            elif np.all([lgp[i_r, i_c]>=180, lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=3600, ts_t10[i_r, i_c]>=3000, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])== True:
                                multi_crop_rain[i_r, i_c] = 3 # Ok
                            
                            elif np.all([lgp[i_r, i_c]>=45, lgp_t5[i_r, i_c]>=120, lgp_t10[i_r, i_c]>=90, ts_t0[i_r, i_c]>=1600, ts_t10[i_r, i_c]>=1200]) == True:
                                multi_crop_rain[i_r, i_c] = 2 #Ok
                            
                            else:
                                multi_crop_rain[i_r, i_c] = 1 #Ok
                            
        
        """Multi cropping zonation for irrigated conditions"""
        for i_r in range(self.im_height):
            for i_c in range(self.im_width):
                
                if self.set_mask:
                    
                    if self.im_mask[i_r, i_c]== self.nodata_val:
                        continue
                    
                    else:
                        
                        if t_climate[i_r, i_c]== 1:
                            
                            if np.all([lgp_t5[i_r, i_c]>=360, lgp_t10[i_r, i_c]>=360, ts_t0[i_r, i_c]>=7200, ts_t10[i_r, i_c]>=7000])==True:
                                multi_crop_irr[i_r, i_c] =8 # ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=300, lgp_t10[i_r, i_c]>=240, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=5100, ts_g_t10[i_r, i_c]>=4800])==True:
                                multi_crop_irr[i_r, i_c] =6 # ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=270, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=5500, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200]) == True:
                                multi_crop_irr[i_r, i_c] =4 # Ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=6400, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_irr[i_r, i_c] =4 #ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])== True:
                                multi_crop_irr[i_r, i_c] =4 # ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=220, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=5500, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700]) == True:
                                multi_crop_irr[i_r, i_c] =3 #Ok
                                
                            elif np.all([lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=6400, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])== True:
                                multi_crop_irr[i_r, i_c] =3 #ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=7200, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])==True:
                                multi_crop_irr[i_r, i_c] =3 # Ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=120, lgp_t10[i_r, i_c]>=90, ts_t0[i_r, i_c]>=1600, ts_t10[i_r, i_c]>=1200]) == True:
                                multi_crop_irr[i_r, i_c] =2 # Ok
                            
                            else:
                                multi_crop_irr[i_r, i_c] =1 # Ok
                        
                        elif t_climate[i_r, i_c] != 1:
                            
                            if np.all([lgp_t5[i_r, i_c]>=360, lgp_t10[i_r, i_c]>=330, ts_t0[i_r, i_c]>=7200, ts_t10[i_r, i_c]>=7000])==True:
                                multi_crop_irr[i_r, i_c] = 8
                            
                            elif np.all([lgp_t5[i_r, i_c]>=330, lgp_t10[i_r, i_c]>=270, ts_t0[i_r, i_c]>=5700, ts_t10[i_r, i_c]>=5500])==True:
                                multi_crop_irr[i_r, i_c] = 7 # ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=300, lgp_t10[i_r, i_c]>=240, ts_t0[i_r, i_c]>=5400, ts_t10[i_r, i_c]>=5100, ts_g_t5[i_r, i_c]>=5100, ts_g_t10[i_r, i_c]>=4800])==True:
                                multi_crop_irr[i_r, i_c] = 6 #ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=270, lgp_t10[i_r, i_c]>=180, ts_t0[i_r, i_c]>=4800, ts_t10[i_r, i_c]>=4500, ts_g_t5[i_r, i_c]>=4300, ts_g_t10[i_r, i_c]>=4000])==True:
                                multi_crop_irr[i_r, i_c] = 5 #ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=240, lgp_t10[i_r, i_c]>=165, ts_t0[i_r, i_c]>=4500, ts_t10[i_r, i_c]>=3600, ts_g_t5[i_r, i_c]>=4000, ts_g_t10[i_r, i_c]>=3200])==True:
                                multi_crop_irr[i_r, i_c] = 4 #ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=200, lgp_t10[i_r, i_c]>=120, ts_t0[i_r, i_c]>=3600, ts_t10[i_r, i_c]>=3000, ts_g_t5[i_r, i_c]>=3200, ts_g_t10[i_r, i_c]>=2700])==True:
                                multi_crop_irr[i_r, i_c] = 3 # ok
                            
                            elif np.all([lgp_t5[i_r, i_c]>=120, lgp_t10[i_r, i_c]>=90, ts_t0[i_r, i_c]>=1600, ts_t10[i_r, i_c]>=1200])==True:
                                multi_crop_irr[i_r, i_c] = 2 #ok
                            
                            else:
                                multi_crop_irr[i_r, i_c] = 1

        if self.set_mask:
            return [np.where(self.im_mask, multi_crop_rain, np.nan), np.where(self.im_mask, multi_crop_irr, np.nan)]
        else:        
            return [multi_crop_rain, multi_crop_irr]
    
    def getAnnualTemperatureAmplitude(self):
        """
        Calculate the annual temperature amplitude (temperature difference from the temperature of the hottest month by temperature from 
        the coldest month).

        Args:
            None.
        Return:
            ann_temp_amp (2D NumPy Array): annual temperature amplitude [Unit: Deg Celsius]
        """

        ann_temp_amp = np.zeros((self.im_height, self.im_width))

        for i in range(self.im_height):
            for j in range(self.im_width):

                if self.set_mask:
                    if self.im_mask[i, j]== self.nodata_val:
                        continue

                monthly_meanT = averageDailyToMonthly(self.meanT_daily[i,j,:], self.leap_year)
                ann_temp_amp[i,j] = np.nanmax(monthly_meanT) - np.nanmin(monthly_meanT)
        
        return ann_temp_amp
    
    def getETODaily(self):
        """
        Get the annual total reference potential evapotranspiration (ET0) calculated from the input climatic variables.
        
        Args:
            None.
        Return:
            ETo (3D NumPy Array): reference evapotranspiration [Unit: mm/day]
        """
        return np.sum(self.pet_daily, axis = 2)
    
    def getAnnualMoistureAvailabilityIndex(self):
        """
        Get the annual moisture availability index (P/ETo * 100). Values ranges from 0 (P < ETo)
        to 100 (P >= ETo). 
        
        Args:
            None.
        Return:
            P/ETO * 100 (2D NumPy Array): annual moisture availability index [Unit: Unitless].
        """
        psum = np.sum(self.totalPrec_daily, axis = 2)
        pet0 = np.sum(self.pet_daily, axis = 2)
        mean_P_ETO = np.zeros(psum.shape)

        mean_P_ETO = np.where(pet0>1e-5, np.round(100.* psum/pet0, 0), np.zeros(pet0.shape))
        mean_P_ETO[mean_P_ETO>1000] = 1000
        mean_P_ETO = np.where(pet0>1e-5, np.round(100.* psum/pet0, 0), np.zeros(pet0.shape))
        mean_P_ETO[np.logical_and(pet0 <1e-5, psum>0)] = 1000.
        mean_P_ETO[np.logical_and(pet0 <1e-5, psum<=0)] = 100.

        # curtailing overshooting values based on Fortran routine  
        mean_P_ETO = np.nanmin([mean_P_ETO, np.full(psum.shape, 30000.)], axis = 0)
        if self.set_mask:
            mean_P_ETO[self.im_mask == self.nodata_val] = 0.
       
        return mean_P_ETO
    
    def getETaDaily(self):
        """
        Get the annual total reference actual evapotranspiration (ETa). ONLY this function can proceed after setting the climatic variables in 
        the object class.
        
        Args:
            None.
        Return:
            eta: (2D-NumPy Array): actual evapotranspiration [Unit: mm].
        """
        return np.sum(self.Eta365, axis = 2)
    
    def getReferenceAnnualWaterDeficit(self):
        """
        Get the annual total reference water deficit (Wde). ONLY this function can proceed after setting the climatic variables in 
        the object class.
        
        Args:a
            None.
        Return:
            ETa (2D-NumPy Array): actual evapotranspiration [Unit: mm].
        """
        eta_sum = np.sum(self.Eta365, axis = 2)
        etm_sum = np.sum(self.Etm365, axis = 2)
        wde = etm_sum - eta_sum

        return wde
    
    def getNPP(self, irr_or_rain:str):
        """
        Calculate the net primary productivity (NPP) for either rainfed or irrigated 
        condition.NPP is estimated as the a functino of incoming solar radiation and 
        soil moisture at the rhizosphere.
        
        Args:
            irr_or_rain (str): Irrigated ('I')/ Rainfed ('R')
        Return:
            NPP (2D NumPy Array): net primary productivity for rainfed/irrigated condition.
        """
        if irr_or_rain == None:
            raise Exception('Please provide string value of I for irrigated or R for rainfed.')

        Rn = np.zeros(366) if self.leap_year else np.zeros(365)
        doy = 366 if self.leap_year else 365
        monthly_shortrad = np.zeros((self.im_height, self.im_width, 12))
        monthly_pr = np.zeros((self.im_height, self.im_width, 12))

        for i in range(self.im_height):
            for j in range(self.im_width):
                Rn = calculateNetRadiationFlux(1, doy, self.latitude[i,j], self.elevation[i,j],  
                                                      self.minT_daily[i,j,:], self.maxT_daily[i,j,:], 
                                                        self.shortrad_daily_MJm2day[i,j,:],
                                                        self.wind_daily[i,j,:], self.rel_humidity_daily[i,j,:],
                                                        self.leap_year)
                
                # Rn = Rn /2.45
                monthly_shortrad[i,j,:] = averageDailyToMonthly(Rn, self.leap_year)
                monthly_pr[i,j,:] = averageDailyToMonthly(self.totalPrec_daily[i,j,:], self.leap_year)


        total_shrad = np.sum(monthly_shortrad, axis =2)
        total_pr = np.sum(monthly_pr, axis = 2)
        
        # Calculate radiative dryness index (Uchijuma and Seino, 1988)
        rdi = np.divide(total_shrad, total_pr, where = total_pr >0, out= np.zeros((self.im_height, self.im_width)))

        if irr_or_rain == 'I':
            # NPP for irrigated condition, ETa = ETm and rdi is 1.375 (maximum of function term)
            rdi = np.nanmin([rdi, np.full((self.im_height, self.im_width), 1.375)], axis = 0)
            npp = np.sum(self.Eto365, axis = 2) * rdi * np.exp(- np.sqrt(9.87+(6.25*rdi))) * 1000
        else:
            # NPP for rainfed condition, ETa is estimated from GAEZ reference water balance
            npp =  np.sum(self.Eta365, axis = 2) * rdi * np.exp(- np.sqrt(9.87+(6.25*rdi))) * 1000
        
        return np.round(npp, 1)

    def getLGPlongest_ver2(self):
        """
        Calculate the total growing days of the longest cycle in a single year.
        
        Args:
            None.
        Return:
            lgp_longest [2-D NumPy Array]: total growing periods of the longest LGP cycle. [Unit: Days]
        """

        lgp_longest = np.zeros((self.im_height, self.im_width), dtype = int)

        eta = self.Eta365.copy()
        etm = self.Etm365.copy()
        Tm = self.meanT_daily.copy()

        for i in range(self.im_height):
            for j in range(self.im_width):

                if self.set_mask:
                    if self.im_mask[i, j]== self.nodata_val:
                        continue
                
                islgp = islgpt(Tm[i,j,:])
                xx = val10day(eta[i,j,:])
                yy = val10day(etm[i,j,:])

                # etamin = np.nanmin(xx)
                # etamax = np.nanmax(xx)

                # etaminidx = np.argmin(xx)
                # etamaxidx = np.argmax(xx)

                # zz = etamin + 

                lgp_whole = np.divide(xx, yy, where= yy>0, out = np.ones(xx.shape))

                count = []

                for k in range(len(lgp_whole)):
                    if islgp[k] == 1 and lgp_whole[k] >=0.4:
                        count.append(1)
                    else:
                        count.append(0)
                

                # find the length of the LGP cycles
                lgp_components = search_cycles(count)
                
                # if there are no growing periods year-round, skip calculation.
                if len(lgp_components[0])==0:
                    lgp_longest[i,j] = 0
                else:
                    # find the longest component
                    sum_list = []
                    
                    for k in range(len(lgp_components[0])):
                        sum_list.append(sum(lgp_components[0][k]))
                    lgp_longest[i,j] = int(np.nanmax(sum_list))
                
        return lgp_longest
    
    def getLGPlongestBeginDate(self):
        """
        Calculate the the beginning day of the longest LGP cycle in a single year frame.
        
        Args:
            None.
        Return:
            lgp_longest_d [2-D NumPy Array]: beginning day of total growing days of the longest. [Unit: DOY]
        """

        lgp_longest_d = np.zeros((self.im_height, self.im_width), dtype = int)

        eta = self.Eta365.copy()
        etm = self.Etm365.copy()
        Tm = self.meanT_daily.copy()

        for i in range(self.im_height):
            for j in range(self.im_width):
                
                # print(f'Row {i}, Col {j}')
                if self.set_mask:
                    if self.im_mask[i, j]== self.nodata_val:
                        continue
                
                islgp = islgpt(Tm[i,j,:])
                xx = val10day(eta[i,j,:])
                yy = val10day(etm[i,j,:])
                lgp_whole = np.divide(xx, yy, where= yy>0, out = np.ones(xx.shape))

                count = []

                for k in range(len(lgp_whole)):
                    if islgp[k] == 1 and lgp_whole[k] >=0.4:
                        count.append(1)
                    else:
                        count.append(0)
                
                # if there are no growing days year-round, skip cycle searching
                if sum(count) ==0:
                    continue
                # find the length of the cycles
                lgp_components = search_cycles(count)

                # if there are no growing periods year-round, skip calculation.
                if len(lgp_components[0])==0:
                    lgp_longest_d[i,j] = 0
                else:
                    # find all days of each cycle
                    sum_list = []
                    
                    for k in lgp_components[0]:
                        if len(k) ==0:
                            sum_list.append(0)
                        else:
                            sum_list.append(sum(k))
                    
                    # change the list into numpy array for possible error occurrence
                    sum_list = np.array(sum_list)
                    lgp_bd = lgp_components[1]
                    idx = np.argwhere(sum_list == np.nanmax(sum_list))[0][0]
                    lgp_longest_d[i,j] = lgp_bd[idx] +1

        return lgp_longest_d
    
    def getBeginningDateofHibernationPeriod(self):
        """
        Calculate the starting date of the hibernation period in a single year frame.
        
        Args:
            None.
        Return:
            lgh_b [2-D NumPy Array]: Beginning date of the hibernation (dormancy period) [Unit: DOY].
        """

        # critical temperature threshold of the least sensitive crop: winter rye
        lgh_b = np.zeros((self.im_height, self.im_width), dtype = int)
        cbtr1, cbtr2 = (-11, -16)

        for i in range(self.im_height):
            for j in range(self.im_width):
                
                # print(f'Row {i}, Col {j}')
                if self.set_mask:
                    if self.im_mask[i, j]== self.nodata_val:
                        continue
                
                # calculate monthly mean average= 
                monthly_Tm = averageDailyToMonthly(self.meanT_daily[i,j,:], self.leap_year)
                tadif0 = np.nanmax(monthly_Tm) - np.nanmin(monthly_Tm)
                
                # Determine the critical breaking temperature (cbtr)
                if tadif0 > 35:
                    cbtr = cbtr1
                elif tadif0 > 20:
                    cbtr = cbtr1 + (cbtr2 - cbtr1) * (35 - tadif0) / 15
                else:
                    cbtr = cbtr2
                
                # three criteria must be satisfied for dormancy (hibernation period) determination
                # 1. Average Temperature must be less than 5
                # 2. Dormancy period must be less than 200 days
                # 3. Average Temperature must satisfy crop-specific temperature threshold.
                dormancy_days = (self.meanT_daily[i,j,:] < 5) & (self.meanT_daily[i,j,:] >= cbtr)

                idx:int = 0
                for k in range(len(dormancy_days)):
                    if dormancy_days[k] == 1:
                        idx = k
                        break
                
                lgh_b[i,j] = idx+1
        
        return lgh_b
    
    def getHibernationPeriodLength(self):
        """
        Calculate the total number of hibernating periods in a single year frame.
        
        Args:
            None.
        Return:
            lgh [2-D NumPy Array]: length of the hibernation period [Unit: Days].
        """

        # critical temperature threshold of the least sensitive crop: winter rye
        lgh = np.zeros((self.im_height, self.im_width), dtype = int)
        cbtr1, cbtr2 = (-11, -16)

        for i in range(self.im_height):
            for j in range(self.im_width):
                
                # print(f'Row {i}, Col {j}')
                if self.set_mask:
                    if self.im_mask[i, j]== self.nodata_val:
                        continue
                
                # calculate monthly mean average= 
                monthly_Tm = averageDailyToMonthly(self.meanT_daily[i,j,:], self.leap_year)
                tadif0 = np.nanmax(monthly_Tm) - np.nanmin(monthly_Tm)
                
                # Determine the critical breaking temperature (cbtr)
                if tadif0 > 35:
                    cbtr = cbtr1
                elif tadif0 > 20:
                    cbtr = cbtr1 + (cbtr2 - cbtr1) * (35 - tadif0) / 15
                else:
                    cbtr = cbtr2
                
                # three criteria must be satisfied for dormancy (hibernation period) determination
                # 1. Average Temperature must be less than 5
                # 2. Dormancy period must be less than 200 days
                # 3. Average Temperature must satisfy crop-specific temperature threshold.
                dormancy_days = (self.meanT_daily[i,j,:] < 5) & (self.meanT_daily[i,j,:] >= cbtr)
                
                # sum up all days with satisfied conditions
                lgh[i,j] = np.nansum(dormancy_days)
        
        return lgh

    ### Dario Spiller additional functions for the revised evaluation of LGP

    def moving_avg_10day(self, x):
        """
        10-day moving average like Fortran val10day(0, MD365, ...)
        - circular convolution over the year (wrap around)
        - returns an array of length MD365
        """
        MD365 = x.shape[-1]
        # duplicate year for wrap-around, then apply window
        x2 = np.concatenate([x, x[:10]], axis=-1)  # extra 10 days for trailing window
        # sliding mean with window size 10, center-aligned to the last day in window
        # emulate Fortran's day i referring to the average ending at i
        out = np.empty(MD365, dtype=np.float32)
        for i in range(MD365):
            # average of days i-9..i (modulo)
            start = i
            end = i + 10
            out[i] = np.nanmean(x2[start:end])
        return out
    
    def normalize_day(self, d, MD365=365):
        """Fortran setdat(): bring day index back into [1..MD365]; Python returns [0..MD365-1]."""
        return ((d % MD365) + MD365) % MD365
    
    def detect_dormancy(self, islgpt):
        """
        Approximate Fortran dormancy detection.
        Returns begdrm (start index), enddrm (end index) in [0..364].
        If no dormancy, returns (None, None).
        """
        MD365 = len(islgpt)
        # We look for a contiguous block of islgpt==0 bounded by 1s
        # Track transitions 1->0 (start) and 0->1 (end)
        prev = islgpt[-1]
        beg = None
        end = None
        for i in range(MD365):
            cur = islgpt[i]
            if prev == 1 and cur == 0 and beg is None:
                beg = i           # first non-growing after growing
            if prev == 0 and cur == 1 and beg is not None and end is None:
                end = i - 1       # last non-growing before growing resumes
            prev = cur
        # if we started dormancy but never ended before year-end, check wrap-around
        if beg is not None and end is None:
            # if year ends in dormancy, end is last 0 before next year's first 1.
            # For simplicity, if entire year is 0s, no dormancy used.
            if np.all(islgpt == 0):
                beg, end = None, None
            else:
                # Find first 1; we already scanned, so if first day is 1 then end is day before beg
                first_one = np.where(islgpt == 1)[0]
                if len(first_one) > 0:
                    end = (beg - 1) % MD365
        return beg, end
    
    def merge_small_gaps(self, components, MD365 = 365, mdbreak = 10):
        """
        Merge consecutive components when the gap between them is <= mdbreak,
        then normalize output to day-of-year (0..MD365-1) with correct (non-negative) lengths.
    
        Parameters
        ----------
        components : list of (beg, end, length, ndwet, ndpet)
            beg, end in 0..MD365-1. If a component spans year end, 'end' < 'beg'.
            'length' in the input will be recomputed from extended indices.
        MD365 : int
            Number of days in the (reference) year. Default 365.
        mdbreak : int
            Maximum gap (days) between components to merge.
    
        Returns
        -------
        merged : list of (beg0, end0, length, ndwet, ndpet)
            beg0/end0 normalized to 0..MD365-1, 'length' computed in extended space (always >= 1).
        """
        if not components:
            return []
    
        # 1) Convert to extended indices (handle wrap-around) and discard input 'length'
        ext = []
        for beg, end, _length_in, ndwet, ndpet in components:
            if end >= beg:
                beg_ext, end_ext = beg, end
            else:
                # spans across year end: extend end by +MD365
                beg_ext, end_ext = beg, end + MD365
            ext.append((beg_ext, end_ext, ndwet, ndpet))
    
        # 2) Sort by extended start and merge small gaps in extended domain
        ext.sort(key=lambda c: c[0])
        merged_ext = [ext[0]]
        for b, e, ndw, ndp in ext[1:]:
            B, E, NDW, NDP = merged_ext[-1]
            gap = b - E - 1
            if gap <= mdbreak:
                # merge with previous: extend end, sum counters
                merged_ext[-1] = (B, e, NDW + ndw, NDP + ndp)
            else:
                merged_ext.append((b, e, ndw, ndp))
    
        # 3) Optional wrap-around merge (last with first) across year end
        if len(merged_ext) > 1:
            B1, E1, NDW1, NDP1 = merged_ext[0]
            BL, EL, NDWL, NDPL = merged_ext[-1]
            gap_wrap = (B1 + MD365) - EL - 1
            if gap_wrap <= mdbreak:
                # merge last->first into one continuous component
                merged_ext = [(BL, E1 + MD365, NDWL + NDW1, NDPL + NDP1)]
    
        # 4) Normalize back to day-of-year and compute length from extended indices
        result = []
        for B, E, NDW, NDP in merged_ext:
            length = (E - B + 1)  # guaranteed non-negative in extended space
            beg0 = B % MD365
            end0 = E % MD365
            result.append((beg0, end0, int(length), int(NDW), int(NDP)))
    
        return result
    
    def getLGPlongest(self,
                      RPlim1=None, RPlim2=None, RPlim3=None,
                      MDBREAK=10, lenmin=30, PHENSTART=0.25):
        """
        Compute the length (days) of the longest growing period (LGP) for each pixel,
        following the IIASA LGP algorithm (Fortran in LGP.txt) once daily balances are available.
    
        Inputs expected on `self`:
          - self.Eta365: (H, W, 365) actual evapotranspiration [mm/day]
          - self.Etm365: (H, W, 365) maximum evapotranspiration [mm/day]
          - self.meanT_daily or self.islgpt: (H, W, 365) temperature-based growing season flag (1/0)
          - Optional: self.totalPrec_daily: (H, W, 365) precipitation [mm/day]
    
        Parameters:
          - RPlim1: start/continue criterion as fraction of ETm for ETa (default deduced or 0.4)
          - RPlim2: end criterion as fraction of ETm for ETa (default deduced or 0.4)
          - RPlim3: rainfall start criterion as fraction of ETm (default deduced or 0.0 if no rainfall)
          - MDBREAK: maximum days between components to be merged (default 10)
          - lenmin: minimum component length to keep (default 30)
          - PHENSTART: 0.25 (phenology start at 25% of ETm range for year-round case)
    
        Returns:
          - lgp_longest: (H, W) int, length (days) of the longest growing period.
          - lgp_beginday: (H, W) int, starting day of the longest growing period.
        """
        H, W, MD365 = self.Eta365.shape
        assert MD365 == 365 or MD365 == 366, "Expected 365-days or 366-days inputs."
    
        # Determine flags: islgpt (1 growing-temperature season, 0 otherwise)
        if hasattr(self, "islgpt365"):
            islgpt_arr = self.islgpt365
        else:
            if not hasattr(self, "meanT_daily"):
                raise ValueError("Provide either self.islgpt365 or self.meanT_daily + islgpt().")
            islgpt_arr = np.zeros_like(self.meanT_daily, dtype=np.int8)
            # Fallback: consider growing if mean T >= 5°C (approximate Fortran Ta >= 5 threshold)
            islgpt_arr = (self.meanT_daily >= 5.0).astype(np.int8)
    
        # Rainfall optional
        has_rain = hasattr(self, "totalPrec_daily") and (self.totalPrec_daily is not None)
    
        # Default thresholds: prefer class/crop-config if present, else safe defaults
        if RPlim1 is None:
            RPlim1 = getattr(self, "RPlim1", 0.4)  # typical crop start/continue
        if RPlim2 is None:
            RPlim2 = getattr(self, "RPlim2", RPlim1)  # end when ETa falls below same fraction
        if RPlim3 is None:
            RPlim3 = getattr(self, "RPlim3", 0.3 if has_rain else 0.0)  # like Fortran: 0 after cold-break
    
        lgp_longest = np.zeros((H, W), dtype=np.int32)
        lgp_beginday = np.zeros((H, W), dtype=np.int32)
    
        # Main loops over pixels (keep clear logic; vectorization is possible later)
        for i in range(H):
            for j in range(W):
                # mask    
                if getattr(self, "set_mask", False):
                    if self.im_mask[i, j] == self.nodata_val:
                        continue
    
                # 10-day averages
                xx = self.moving_avg_10day(self.Eta365[i, j, :])  # ETa 10-day
                yy = self.moving_avg_10day(self.Etm365[i, j, :])  # ETm 10-day
                if has_rain:
                    zz = self.moving_avg_10day(self.totalPrec_daily[i, j, :])  # rain 10-day
                else:
                    zz = np.zeros_like(xx)
    
                islgp = islgpt_arr[i, j, :].astype(np.int8)
    
                # Duplicate the year for scanning across boundaries
                xx2 = np.concatenate([xx, xx])       # length 730
                yy2 = np.concatenate([yy, yy])
                zz2 = np.concatenate([zz, zz])
                islgp2 = np.concatenate([islgp, islgp])
    
                # ---- Fortran: scan for first break day to set scanning window (istrt0, istrt1) ----
                # "break" means a day failing season or ETa < RPlim2 * ETm
                istrt0 = None
                for d in range(MD365):
                    if (islgp[d] == 0) or (yy[d] > 0 and xx[d] < yy[d] * RPlim2):
                        istrt0 = (d + 1) % MD365
                        break
                if istrt0 is None:
                    # Year-round case: no break
                    # Fortran: set lgp = full year and pick phenology start at 25% of ETm range
                    ETmin = np.nanmin(yy)
                    ETmax = np.nanmax(yy)
                    zz_thr = ETmin + PHENSTART * (ETmax - ETmin)
    
                    # Find first day after ETm rises above threshold scanning from min position
                    jETmn = int(np.nanargmin(yy))
                    beglgp = None
                    for k in range(jETmn, jETmn + MD365):
                        if yy[k % MD365] >= zz_thr:
                            beglgp = k % MD365
                            break
                    # If humid all-year you may set LGP=366 in Fortran; here we keep MD365.
                    lgp_longest[i, j] = MD365
                    continue
    
                # scanning window over two-year sequence
                start_idx = istrt0      # 0..364
                end_idx = start_idx + MD365 - 1
                
                # ---- Fortran: dormancy detection and rainfall-start suppression right after dormancy ----
                begdrm, enddrm = self.detect_dormancy(islgp)
                # map to the 730-day window indexes if present
                enddrm0 = None if enddrm is None else enddrm
                components = []
                in_lgp = False
                cur_beg = None
                ndwet = 0   # rainy >= ETm
                ndpet = 0   # ETa >= ETm
    
                # Iterate over 2-year window (indices in 0..729)
                for t in range(start_idx, end_idx):  # MD365 days scanned within 2-year arrays - np.minimum(end_idx + MD365, MD365*2 - 1)
                    ii = t  # absolute index
                    day_mod = ii % MD365
    
                    xx_t = xx2[ii]
                    yy_t = yy2[ii]
                    zz_t = zz2[ii]
                    is_season = islgp2[ii] == 1
    
                    # Determine rainfall start limit (RPl3): suppress right after dormancy end
                    if begdrm is not None and enddrm is not None:
                        after_cold_break = (day_mod == (enddrm + 1) % MD365) or (day_mod == (enddrm0 + 1) % MD365)
                    else:
                        after_cold_break = False
                    RPl3_use = 0.0 if after_cold_break else RPlim3
    
                    # Check LGP start condition (Fortran: islgpt==1 and ETa>=RPlim1*ETm and rain>=RPl3*ETm)
                    start_ok = (is_season and
                                (yy_t <= 0 or xx_t >= yy_t * RPlim1) and
                                (yy_t <= 0 or zz_t >= yy_t * RPl3_use))
    
                    end_ok = (not is_season) or (yy_t > 0 and xx_t < yy_t * RPlim2)
    
                    if not in_lgp and start_ok:
                        in_lgp = True
                        cur_beg = ii
                        ndwet = 0
                        ndpet = 0
    
                    elif in_lgp and end_ok:
                        in_lgp = False
                        cur_end = ii - 1
                        # record component in day-of-year coordinates (0..364)
                        beg_d = cur_beg % MD365
                        end_d = cur_end % MD365
                        # length across possible wrap inside window
                        length = (cur_end - cur_beg + 1)
                        components.append((beg_d, end_d, length, ndwet, ndpet))
    
                    # accumulate day types while inside LGP
                    if in_lgp:
                        if yy_t > 0 and zz_t >= yy_t:   # rainy day >= ETm
                            ndwet += 1
                            ndpet += 1
                        elif yy_t > 0 and xx_t >= yy_t: # ETa >= ETm (PET-satisfied)
                            ndpet += 1
    
                # Close trailing component if we finished inside LGP
                if in_lgp and cur_beg is not None:
                    cur_end = end_idx + MD365 - 1
                    beg_d = cur_beg % MD365
                    end_d = cur_end % MD365
                    length = (cur_end - cur_beg + 1) % MD365
                    components.append((beg_d, end_d, length, ndwet, ndpet))

                # ---- Fortran: merge components with small gaps (MDBREAK) including wrap-around ----
                components = self.merge_small_gaps(components, MD365=MD365, mdbreak=MDBREAK)

                # ---- Fortran: discard components shorter than lenmin ----
                components = [c for c in components if c[2] >= lenmin]

                if not components:
                    lgp_longest[i, j] = 0
                    continue
    
                # ---- Fortran: sort by length (descending) and pick longest ----
                components.sort(key=lambda c: c[2], reverse=True)
                longest = components[0]
                lgp_longest[i, j] = int(longest[2])
                lgp_beginday[i, j] = int(longest[0])
                
        return lgp_longest, lgp_beginday


#----------------- End of file -------------------------#


