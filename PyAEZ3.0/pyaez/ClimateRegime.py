"""
PyAEZ Version 4.0 (March 2026)

The `ClimateRegime` class is responsible for reading, loading, and calculating 
the agro-climatic indicators required to run PyAEZ.

Authors and Contributors:
- 2021: N. Lakmal Deshapriya
- 2022–2023: Swun Wunna Htet and Kittiphon Boonma
- 2024 (April): Swun Wunna Htet (up to version 2.3)
- 2025 (October): Dario Spiller
- 2026 (March): Riley Tuccio

Modification History:

Up to Version 2.3:
1. Removed object class declarations from other modules.
2. Integrated reference water balance calculations into the routine.
3. Added new agro-climatic indicator functions.
4. Enhanced handling of monthly vs. daily time dimensions.

From Version 3.0:
1. Improved leap year management.
2. Reviewed and refined water balance calculations.

Version 4.0: 
1. Implemented Vectorization for raster computations
2. Implemented Chunking Calculations for best performance for functions
3. consalidated conditional code of if-elif-else statements into condition checks for masking
"""

import numpy as np

from UtilitiesCalc import compute_chunk_size_multi, generateLatitudeMap, interpMonthlyToDaily, averageDailyToMonthly
from ETOCalc import calculateETONumba, compute_Rn_chunk
from LGPCalc import psh, RefWaterBalanceCalc, rainPeak, islgpt, val10day, search_cycles, process_chunk
from ThermalScreening import getTempTrend, getSmoothTemp, getTemperatureGrowingPeriod
from typing import List, Tuple
import psutil


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


        # flattening data for easier vectorization
        n_pixels = self.im_height * self.im_width

        #Initialize daily data arrays
        self.meanT_daily = np.zeros((n_pixels, doy))
        self.totalPrec_daily = np.zeros((n_pixels, doy))
        self.minT_daily = np.zeros((n_pixels, doy))
        self.maxT_daily = np.zeros((n_pixels, doy))
        self.shortrad_daily = np.zeros((n_pixels, doy))
        self.wind_daily = np.zeros((n_pixels, doy))
        self.rel_humidity_daily = np.zeros((n_pixels, doy))
        self.pet_daily = np.zeros((n_pixels, doy))
        self.shortrad_daily_MJm2day = np.zeros((n_pixels, doy))

        #flattening inputs
        minT = min_temp.reshape(n_pixels, -1)
        maxT = max_temp.reshape(n_pixels, -1)
        prec = precipitation.reshape(n_pixels, -1)
        short = short_rad.reshape(n_pixels, -1)
        wind = wind_speed.reshape(n_pixels, -1)
        humid = rel_humidity.reshape(n_pixels, -1)
        lat = self.latitude.reshape(n_pixels)
        elev = self.elevation.reshape(n_pixels)



        ## Calcuating chunk size to be used

        chunk_size = compute_chunk_size_multi(
            arrays=[minT, maxT, prec, short, wind, humid, lat, elev],
            extra_arrays= None)

        if self.set_mask:
            mask = (self.im_mask != self.nodata_val).reshape(n_pixels)
        else:
            mask = np.ones(n_pixels, dtype=bool)

        # --- monthly mean ---
        meanT_monthly = (minT + maxT) / 2

        # --- chunk loop ---
        for start in range(0, n_pixels, chunk_size):
            end = min(start + chunk_size, n_pixels)

            chunk_mask = mask[start:end]
            if not np.any(chunk_mask):
                continue

            idx = np.where(chunk_mask)[0] + start  # global indices

            #Instead of doing looping by i_row, i_col, this code does it by the pixels within the chunk of memory calculated
            if self.set_monthly:
                meanT_chunk = interpMonthlyToDaily(meanT_monthly[idx], 1, doy)
                prec_chunk = interpMonthlyToDaily(prec[idx], 1, doy, no_minus_values=True)
                minT_chunk = interpMonthlyToDaily(minT[idx], 1, doy)
                maxT_chunk = interpMonthlyToDaily(maxT[idx], 1, doy)
                short_chunk = interpMonthlyToDaily(short[idx], 1, doy, no_minus_values=True)
                wind_chunk = interpMonthlyToDaily(wind[idx], 1, doy, no_minus_values=True)
                humid_chunk = interpMonthlyToDaily(humid[idx], 1, doy, no_minus_values=True)
            else:
                meanT_chunk = (minT[idx] + maxT[idx]) / 2
                prec_chunk = prec[idx]
                minT_chunk = minT[idx]
                maxT_chunk = maxT[idx]
                short_chunk = short[idx]
                wind_chunk = wind[idx]
                humid_chunk = humid[idx]

            # radiation conversion (vectorized)
            short_MJ_chunk = short_chunk * 86400.0 / 1e6


            # STORE (vectorized)
            self.meanT_daily[idx] = meanT_chunk
            self.totalPrec_daily[idx] = prec_chunk
            self.minT_daily[idx] = minT_chunk
            self.maxT_daily[idx] = maxT_chunk
            self.shortrad_daily[idx] = short_chunk
            self.wind_daily[idx] = wind_chunk
            self.rel_humidity_daily[idx] = humid_chunk
            self.shortrad_daily_MJm2day[idx] = short_MJ_chunk



            for k, p in enumerate(idx):
                self.pet_daily[p] = calculateETONumba(
                    1, doy,
                    lat[p],
                    elev[p],
                    minT_chunk[k],
                    maxT_chunk[k],
                    wind_chunk[k],
                    short_MJ_chunk[k],
                    humid_chunk[k],
                    self.leap_year
                )

        # --- reshape back to (H, W, time) ---
        self.meanT_daily = self.meanT_daily.reshape(self.im_height, self.im_width, doy)
        self.totalPrec_daily = self.totalPrec_daily.reshape(self.im_height, self.im_width, doy)
        self.minT_daily = self.minT_daily.reshape(self.im_height, self.im_width, doy)
        self.maxT_daily = self.maxT_daily.reshape(self.im_height, self.im_width, doy)
        self.shortrad_daily = self.shortrad_daily.reshape(self.im_height, self.im_width, doy)
        self.wind_daily = self.wind_daily.reshape(self.im_height, self.im_width, doy)
        self.rel_humidity_daily = self.rel_humidity_daily.reshape(self.im_height, self.im_width, doy)
        self.pet_daily = self.pet_daily.reshape(self.im_height, self.im_width, doy)
        self.shortrad_daily_MJm2day = self.shortrad_daily_MJm2day.reshape(self.im_height, self.im_width, doy)

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
        Txsnm = 0.0  # Snow melt temperature threshold (°C)
        Fsnm = 5.5  # Snow melting coefficient

        # Variables initialization - Aliases (reduce attribute lookups)
        Tx365 = self.maxT_daily  # shape: (H, W, D)
        Ta365 = self.meanT_daily
        Pcp365 = self.totalPrec_daily
        self.Eto365 = self.pet_daily  # Eto
        Eto365 = self.Eto365

        # --- Shapes & pre-allocation ---
        H, W, T = Tx365.shape
        self.Etm365 = np.zeros((H, W, T), dtype=np.float64)
        self.Eta365 = np.zeros((H, W, T), dtype=np.float64)
        self.Sb365 = np.zeros((H, W, T), dtype=np.float64)
        self.Wb365 = np.zeros((H, W, T), dtype=np.float64)
        self.Wx365 = np.zeros((H, W, T), dtype=np.float64)
        self.kc365 = np.zeros((H, W, T), dtype=np.float64)

        # Build mask
        if self.set_mask:
            mask = (self.im_mask != self.nodata_val)
        else:
            mask = np.ones((H, W), dtype=bool)

        rows, cols = np.where(mask)
        n_pix = rows.size

        T = Tx365.shape[2]

        # Flattening Inputs
        Tx_all = Tx365[rows, cols, :].astype(np.float64)
        Ta_all = Ta365[rows, cols, :].astype(np.float64)
        P_all = Pcp365[rows, cols, :].astype(np.float64)
        Eto_all = Eto365[rows, cols, :].astype(np.float64)

        # Growing season length
        lgpt5_all = np.count_nonzero(Ta_all >= 5.0, axis=1).astype(np.int64)

        # Allocate
        istart0_all = np.empty(n_pix, dtype=np.int64)
        istart1_all = np.empty(n_pix, dtype=np.int64)
        istup_all = np.empty((n_pix, T), dtype=np.int64)

        # Compute per-pixel before vectorization with chunking
        for i in range(n_pix):

            Ta_p = Ta_all[i]

            # rainPeak → start/end of growing season
            s0, s1 = rainPeak(Ta_p, lgpt5_all[i])
            istart0_all[i] = s0
            istart1_all[i] = s1

            # temperature trend
            trend = getTempTrend(Ta_p)

            if trend is None or len(trend) < T:
                raise ValueError(f"Invalid istup at pixel {i}")

            istup_all[i, :] = np.asarray(trend, dtype=np.int64)


        # Calculating Optimized Chunk Size
        chunk_size = compute_chunk_size_multi(
            arrays=[Tx_all, Ta_all, P_all, Eto_all, istart0_all, istart1_all, istup_all],
            extra_arrays=None)

        for start in range(0, n_pix, chunk_size):
            end = min(start + chunk_size, n_pix)

            # Slice chunk
            Tx = Tx_all[start:end]
            Ta = Ta_all[start:end]
            P = P_all[start:end]
            Eto = Eto_all[start:end]

            istart0 = istart0_all[start:end]
            istart1 = istart1_all[start:end]
            istup = istup_all[start:end]
            lgpt5 = lgpt5_all[start:end]

            n = Tx.shape[0]

            # arrays to story chunk outputs
            Eta = np.empty((n, T), dtype=np.float64)
            Etm = np.empty((n, T), dtype=np.float64)
            Wb = np.empty((n, T), dtype=np.float64)
            Wx = np.empty((n, T), dtype=np.float64)
            Sb = np.empty((n, T), dtype=np.float64)
            kc = np.empty((n, T), dtype=np.float64)


            ''' this calls to the newly created process_chunk @jit functino in LGP, which handles variable passing 
                in a way that doesn't give nonpython type interference errors, which is a common thing with numba functions.
            '''
            process_chunk(
                Tx, Ta, P, Eto,
                istart0, istart1, istup, lgpt5,
                Txsnm, Fsnm, Sa, D, itflg, T,
                Eta, Etm, Wb, Wx, Sb, kc
            )

           # writing back to the full grid
            r = rows[start:end]
            c = cols[start:end]

            self.Eta365[r, c, :] = Eta
            self.Etm365[r, c, :] = Etm
            self.Wb365[r, c, :] = Wb
            self.Wx365[r, c, :] = Wx
            self.Sb365[r, c, :] = Sb
            self.kc365[r, c, :] = kc



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
            - 
        """
        H, W = self.im_height, self.im_width #Height and Width of the raster
        thermal_climate = np.zeros((H, W), dtype=np.int8)

        chunk_size = compute_chunk_size_multi(
            arrays=[thermal_climate],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # meanT_sl_m
                {"shape_factor": (12,), "dtype": np.float64},  # meanT_m
                {"shape_factor": (12,), "dtype": np.float64},  # P_PET_m
                {"shape_factor": (12,), "dtype": np.float64},  # prec_m

                # Reductions
                {"shape_factor": (1,), "dtype": np.float64},  # Ta_diff
                {"shape_factor": (1,), "dtype": np.float64},  # minT_sl
                {"shape_factor": (1,), "dtype": np.float64},  # meanT_avg
                {"shape_factor": (1,), "dtype": np.int64},  # months_ge10
                {"shape_factor": (1,), "dtype": np.float64},  # total_precip
                {"shape_factor": (1,), "dtype": np.float64},  # rainfall_range

                # Seasonal accumulators
                {"shape_factor": (1,), "dtype": np.float64},  # summer_PET0
                {"shape_factor": (1,), "dtype": np.float64},  # winter_PET0

                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )




        if self.set_mask:
            valid_mask = self.im_mask != self.nodata_val
        else:
            valid_mask = np.ones((H, W), dtype=np.bool_)

        for i0 in range(0, H, chunk_size): #chunk index for height
            for j0 in range(0, W, chunk_size): #chunk index for width
                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)


                mask_chunk = valid_mask[i0:i1, j0:j1]
                if not np.any(mask_chunk):
                    continue

                # Slice data over the chunk values for daily
                lat_c = self.latitude[i0:i1, j0:j1] 
                meanT_sl_c = self.meanT_daily_sealevel[i0:i1, j0:j1, :] 
                meanT_c = self.meanT_daily[i0:i1, j0:j1, :]
                P_PET_c = self.P_by_PET_daily[i0:i1, j0:j1, :]
                prec_c = self.totalPrec_daily[i0:i1, j0:j1, :]

                # Allocate monthly arrays
                hC, wC = i1 - i0, j1 - j0

                meanT_sl_m = np.zeros((hC, wC, 12))
                meanT_m = np.zeros((hC, wC, 12))
                P_PET_m = np.zeros((hC, wC, 12))
                prec_m = np.zeros((hC, wC, 12))

                meanT_sl_m = averageDailyToMonthly(meanT_sl_c, self.leap_year)
                meanT_m = averageDailyToMonthly(meanT_c, self.leap_year)
                P_PET_m = averageDailyToMonthly(P_PET_c, self.leap_year)
                prec_m = averageDailyToMonthly(prec_c, self.leap_year)


                # ignoring non mask valid chunks
                meanT_sl_m[~mask_chunk] = 0
                meanT_m[~mask_chunk] = 0
                P_PET_m[~mask_chunk] = 0
                prec_m[~mask_chunk] = 0


                #Vectorized computations
                Ta_diff = np.max(meanT_sl_m, axis=2) - np.min(meanT_sl_m, axis=2)
                minT_sl = np.min(meanT_sl_m, axis=2)
                meanT_avg = np.mean(meanT_m, axis=2)
                months_ge10 = np.sum(meanT_sl_m >= 10.0, axis=2)

                total_precip = np.sum(prec_c, axis=2)

                rainfall_range = np.max(prec_m, axis=2) - np.min(prec_m, axis=2)

                # Hemisphere masks
                north = lat_c > 5
                south = lat_c < -5
                equatorial = np.abs(lat_c) <= 5

                # Precompute seasonal sums
                summer_idx_n = np.array([3, 4, 5, 6, 7, 8])
                winter_idx_n = np.array([9, 10, 11, 0, 1, 2])
                summer_idx_s = winter_idx_n
                winter_idx_s = summer_idx_n

                summer_PET0 = np.zeros((hC, wC))
                winter_PET0 = np.zeros((hC, wC))


                summer_PET0[north] = np.sum(P_PET_m[north][:, summer_idx_n], axis=1)
                winter_PET0[north] = np.sum(P_PET_m[north][:, winter_idx_n], axis=1)

                summer_PET0[south] = np.sum(P_PET_m[south][:, summer_idx_s], axis=1)
                winter_PET0[south] = np.sum(P_PET_m[south][:, winter_idx_s], axis=1)

                # --- Classification (still vectorized masks where possible) ---
                out = np.zeros((hC, wC), dtype=np.int8)

                # Equatorial
                eq_uniform = equatorial & (rainfall_range < 50)
                eq_seasonal = equatorial & (rainfall_range >= 50)

                out[eq_uniform & (out == 0)] = 1
                out[eq_seasonal & (out == 0)] = 3




                # Tropical
                trop = (minT_sl >= 18.) & (Ta_diff < 15.) & (~equatorial)
                out[trop & (meanT_avg < 20.) & (out==0)] = 2
                out[trop & (meanT_avg >= 20.) & (out==0)] = 1

                # Subtropical
                subtrop = (
                        (minT_sl >= 5.) &
                        (months_ge10 >= 8) &
                        (~trop) &
                        (~equatorial)
                )

                low_prec = subtrop & (total_precip < 250)
                out[low_prec & (out==0)] = 5

                summer_dom = summer_PET0 >= winter_PET0

                out[subtrop & (~low_prec) & summer_dom & (out==0)] = 3
                out[subtrop & (~low_prec) & (~summer_dom) & (out ==0)] = 4

                # Temperate
                temp = (
                        (months_ge10 >= 4) &
                        (~subtrop) &
                        (~trop) &
                        (~equatorial)
                )

                out[temp & (Ta_diff <= 20) & (out==0)] = 6
                out[temp & (Ta_diff > 20) & (Ta_diff <= 35) & (out==0)] = 7
                out[temp & (Ta_diff > 35) & (out==0)] = 8

                # Boreal
                boreal = (
                        (months_ge10 >= 1) &
                        (months_ge10 < 4)
                )

                out[boreal & (Ta_diff <= 20) & (out==0)] = 9
                out[boreal & (Ta_diff > 20) & (Ta_diff <= 35) & (out==0)] = 10
                out[boreal & (Ta_diff > 35) & (out==0)] = 11

                # Arctic
                arctic = months_ge10 == 0
                out[arctic & (out==0)] = 12

                # Apply mask + write back
                thermal_climate[i0:i1, j0:j1] = out * mask_chunk

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

        H, W = self.im_height, self.im_width

        thermal_zone = np.zeros((H, W), dtype=np.float64)

        # mask
        if self.set_mask:
            valid_mask = self.im_mask != self.nodata_val
        else:
            valid_mask = np.ones((H, W), dtype=np.bool_)

        chunk_size = compute_chunk_size_multi(
            arrays=[thermal_zone],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # meanT_c
                {"shape_factor": (12,), "dtype": np.float64},  # meanT_sl_m

                # Reductions
                {"shape_factor": (1,), "dtype": np.float64},  # minT_sl
                {"shape_factor": (1,), "dtype": np.float64},  # maxT
                {"shape_factor": (1,), "dtype": np.float64},  # minT
                {"shape_factor": (1,), "dtype": np.int64},  # months_gt10_sl
                {"shape_factor": (1,), "dtype": np.float64},  # months_lt10

                # out
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )

        for i0 in range(0, H, chunk_size):
            for j0 in range(0, W, chunk_size):
                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)

                mask_chunk = valid_mask[i0:i1, j0:j1]
                if not np.any(mask_chunk):
                    continue

                meanT_c = self.meanT_daily[i0:i1, j0:j1, :]
                meanT_sl_c = self.meanT_daily_sealevel[i0:i1, j0:j1, :]

                hC, wC = i1 - i0, j1 - j0


                meanT_m = np.zeros((hC, wC, 12))
                meanT_sl_m = np.zeros((hC, wC, 12))

                meanT_m  = averageDailyToMonthly(meanT_c,self.leap_year)
                meanT_sl_m = averageDailyToMonthly(meanT_sl_c,self.leap_year)

                meanT_m[~mask_chunk] = 0
                meanT_sl_m[~mask_chunk] = 0

                # --- vectorized features ---
                minT_sl = np.min(meanT_sl_m, axis=2)
                maxT = np.max(meanT_m, axis=2)
                minT = np.min(meanT_m, axis=2)
                meanT_avg = np.mean(meanT_m, axis=2)

                months_gt10_sl = np.sum(meanT_sl_m > 10, axis=2)
                months_ge10_sl = np.sum(meanT_sl_m >= 10, axis=2)

                months_lt5 = np.sum(meanT_m < 5, axis=2)
                months_gt10 = np.sum(meanT_m > 10, axis=2)
                months_lt10 = np.sum(meanT_m < 10, axis=2)

                Ta_diff = maxT - minT

                out = np.zeros((hC, wC), dtype=np.float64)

                # --- Tropics ---
                trop = (minT_sl >= 18) & (Ta_diff < 15)
                out[trop & (meanT_avg > 20)] = 1
                out[trop & (meanT_avg <= 20)] = 2

                # --- Subtropics ---
                subtrop = (minT_sl > 5) & (months_gt10_sl >= 8) & (~trop)

                out[subtrop & (months_lt5 >= 1) & (months_gt10 >= 4)] = 4
                out[subtrop & (months_lt5 >= 1) & (months_gt10 >= 1) &
                    ~(months_gt10 >= 4)] = 5
                out[subtrop & (months_lt10 == 12)] = 6
                out[subtrop & (out == 0)] = 3

                # --- Temperate ---
                temp = (months_ge10_sl >= 4) & (~subtrop) & (~trop)

                out[temp & (months_lt5 >= 1) & (months_gt10 >= 4)] = 7
                out[temp & (months_lt5 >= 1) & (months_gt10 >= 1) &
                    ~(months_gt10 >= 4)] = 8
                out[temp & (months_lt10 == 12)] = 9

                # --- Boreal ---
                boreal = (months_ge10_sl >= 1) & (months_ge10_sl < 4)

                out[boreal & (months_lt5 >= 1) & (months_gt10 >= 1)] = 10
                out[boreal & (months_lt10 == 12)] = 11

                # --- Arctic ---
                arctic = months_ge10_sl == 0
                out[arctic] = 12

                # write back to final output
                thermal_zone[i0:i1, j0:j1] = out * mask_chunk
    
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

        H, W, T = self.meanT_daily.shape
        DAYS_IN_YEAR = 366 if self.leap_year else 365
        days = np.arange(DAYS_IN_YEAR)


        interp_daily_temp = np.zeros((H, W, DAYS_IN_YEAR), dtype=np.float32)

        chunk_size = compute_chunk_size_multi(
            arrays=[interp_daily_temp],
            extra_arrays=None)

        for i0 in range(0, H, chunk_size):
            for j0 in range(0, W, chunk_size):
                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)

                temp_chunk = self.meanT_daily[i0:i1, j0:j1, :]  # (hC, wC, T)
                hC, wC = temp_chunk.shape[:2]

                # reshape → (N_pixels, T)
                temp_2d = temp_chunk.reshape(-1, T)

                # vectorized polynomial fit
                coeffs = np.polyfit(days, temp_2d.T, 5)
                # shape: (deg+1, N_pixels)

                # evaluate polynomial
                interp_2d = np.polyval(coeffs, days[:, None]).T  # (N_pixels, T)

                # reshape back
                interp_chunk = interp_2d.reshape(hC, wC, T)

                interp_daily_temp[i0:i1, j0:j1, :] = interp_chunk


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

        H, W, T = self.Etm365.shape

        lgp_tot = np.zeros((H, W), dtype=np.float32)

        if self.set_mask:
            valid_mask = self.im_mask != self.nodata_val
        else:
            valid_mask = np.ones((H, W), dtype=np.bool_)

        #computing optimized chunk size
        chunk_size = compute_chunk_size_multi(
            arrays=[lgp_tot],
            extra_arrays=None)


        for i0 in range(0, H, chunk_size):
            for j0 in range(0, W, chunk_size):
                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)

                mask_chunk = valid_mask[i0:i1, j0:j1]
                if not np.any(mask_chunk):
                    continue

                Etm_c = self.Etm365[i0:i1, j0:j1, :]
                Eta_c = self.Eta365[i0:i1, j0:j1, :]
                T_c = self.meanT_daily[i0:i1, j0:j1, :]

                hC, wC = Etm_c.shape[:2]

                # reshape → (N_pixels, T)
                Etm_2d = Etm_c.reshape(-1, T)
                Eta_2d = Eta_c.reshape(-1, T)
                T_2d = T_c.reshape(-1, T)
                mask_1d = mask_chunk.reshape(-1)

                # filter valid pixels only
                Etm_2d = Etm_2d[mask_1d]
                Eta_2d = Eta_2d[mask_1d]
                T_2d = T_2d[mask_1d]



                # still per-row unless vectorized
                n_pix = Etm_2d.shape[0]

                xx = np.zeros((n_pix, T), dtype=np.float32)
                yy = np.zeros((n_pix, T), dtype=np.float32)
                islgp_arr = np.zeros((n_pix, T), dtype=np.int8)

                for k in range(n_pix):
                    xx[k, :] = val10day(Eta_2d[k])
                    yy[k, :] = val10day(Etm_2d[k])
                    islgp_arr[k, :] = islgpt(T_2d[k])

                # compute ratio
                lgp_ratio = xx/ yy

                # vectorized condition
                cond = (islgp_arr == 1) & (lgp_ratio >= 0.4)

                counts = np.sum(cond, axis=1)

                # write back
                out_chunk = np.zeros(hC * wC, dtype=np.float32)
                out_chunk[mask_1d] = counts

                lgp_tot[i0:i1, j0:j1] = out_chunk.reshape(hC, wC)


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
        Compute thermal zones for fallow requirements.

        Args:
            tzone (2D np.array): thermal zone classes (from TZone classification)

        Returns:
            2D np.array: tzone_fallow (same shape as input raster)
        """

        # the algorithm needs to calculate the annual mean temperature.

        H, W = tzone.shape
        out = np.zeros((H, W), dtype=int)

        chunk_size = compute_chunk_size_multi(
            arrays=[out],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # annual_Tmean
                {"shape_factor": (12,), "dtype": np.float64},  # monthly_means

                # Seasonal accumulators
                {"shape_factor": (1,), "dtype": np.float64},  # warmest_month

                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],)


        for i in range(0, H, chunk_size):
            for j in range(0, W, chunk_size):
                i_end = min(i + chunk_size, H)
                j_end = min(j + chunk_size, W)

                # Slice chunk
                tzone_chunk = tzone[i:i_end, j:j_end]
                temp_chunk = self.meanT_daily[i:i_end, j:j_end, :]

                #Vectorized computation on chunk
                annual_Tmean = np.mean(temp_chunk, axis=2)
                monthly_means = averageDailyToMonthly(temp_chunk, self.leap_year)
                warmest_month = np.max(monthly_means, axis=2)

                tzonefallow_chunk = np.zeros_like(tzone_chunk, dtype=int)

                tropics = (tzone_chunk == 1) | (tzone_chunk == 2)
                non_tropics = ~tropics

                # Tropical
                tzonefallow_chunk[tropics & (annual_Tmean > 25)] = 1
                tzonefallow_chunk[tropics & (annual_Tmean > 20) & (annual_Tmean <= 25)] = 2
                tzonefallow_chunk[tropics & (annual_Tmean > 15) & (annual_Tmean <= 20)] = 3
                tzonefallow_chunk[tropics & (annual_Tmean <= 15)] = 4

                # Non-tropical
                tzonefallow_chunk[non_tropics & (warmest_month > 20)] = 5
                tzonefallow_chunk[non_tropics & (warmest_month <= 20)] = 6

                # write back to total output 'out'
                out[i:i_end, j:j_end] = tzonefallow_chunk

        if self.set_mask:
            out = np.where(self.im_mask != self.nodata_val, out, np.nan)

        return out


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

        H, W, T = self.meanT_daily.shape

        fi = np.zeros((H, W), dtype=np.float32)
        permafrost = np.zeros((H, W), dtype=np.int8)

        # mask
        if self.set_mask:
            valid_mask = self.im_mask != self.nodata_val
        else:
            valid_mask = np.ones((H, W), dtype=np.bool_)


        chunk_size = compute_chunk_size_multi(
            arrays=[fi,permafrost],
            extra_arrays=None)


        for i0 in range(0, H, chunk_size):
            for j0 in range(0, W, chunk_size):
                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)

                mask_chunk = valid_mask[i0:i1, j0:j1]
                if not np.any(mask_chunk):
                    continue

                temp_c = self.meanT_daily[i0:i1, j0:j1, :]

                # vectorized thawing/freezing indices
                meanT_gt_0 = np.maximum(temp_c, 0)
                meanT_le_0 = np.minimum(temp_c, 0)

                ddt = np.sum(meanT_gt_0, axis=2)
                ddf = -np.sum(meanT_le_0, axis=2)

                # avoid division issues safely
                sqrt_ddt = np.sqrt(ddt)
                sqrt_ddf = np.sqrt(ddf)

                fi_chunk = sqrt_ddf / (sqrt_ddf + sqrt_ddt)

                # replace NaNs
                fi_chunk = np.nan_to_num(fi_chunk)

                # fully vectorized classification
                pf_chunk = np.zeros_like(fi_chunk, dtype=np.int8)

                pf_chunk[fi_chunk > 0.625] = 1
                pf_chunk[(fi_chunk > 0.57) & (fi_chunk <= 0.625)] = 2
                pf_chunk[(fi_chunk > 0.495) & (fi_chunk <= 0.57)] = 3
                pf_chunk[fi_chunk <= 0.495] = 4

                # apply mask
                fi[i0:i1, j0:j1] = fi_chunk * mask_chunk
                permafrost[i0:i1, j0:j1] = pf_chunk * mask_chunk

        if self.set_mask:
            return [
                np.where(self.im_mask, fi, np.nan),
                np.where(self.im_mask, permafrost, np.nan)
            ]
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
        
        # reclassifying the existing 12 classes of thermal climate into 6 major thermal climate.
        # Class 1: Tropics, lowland
        # Class 2: Tropics, highland
        # Class 3: Subtropics
        # Class 4: Temperate Climate
        # Class 5: Boreal Climate
        # Class 6: Arctic Climate

        nodata_val = self.nodata_val
        H, W = self.im_height, self.im_width

        # Reclassing tclimate
        lut = np.zeros(13, dtype=np.int8)
        lut[1] = 1
        lut[2] = 2
        lut[3] = 3
        lut[4] = 3
        lut[5] = 3
        lut[6:9] = 4 # 6,7,8 → 4
        lut[9:12] = 5 # 9,10,11 → 5
        lut[12] = 6

        aez_tclimate = np.zeros((H, W), dtype=np.int8)
        aez_temp_regime = np.zeros((H, W), dtype=np.int8)
        aez_moisture_regime = np.zeros((H, W), dtype=np.int8)
        aez = np.zeros((H, W), dtype=np.int8)

        chunk_size = compute_chunk_size_multi(
            arrays=[aez_tclimate, aez_temp_regime, aez_moisture_regime, aez],
            extra_arrays=None
        )

        # ---- Step 1a: Thermal climate ----
        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                t_chunk = tclimate[r0:r1, c0:c1]
                t_chunk_int = np.where(np.isnan(t_chunk), 0, t_chunk)
                t_chunk_int = np.clip(t_chunk_int.astype(np.int8), 0, len(lut) - 1)

                out_chunk = lut[t_chunk_int]

                if self.set_mask:
                    mask_chunk = self.im_mask[r0:r1, c0:c1] != nodata_val
                    out_chunk = np.where(mask_chunk, out_chunk, nodata_val)

                aez_tclimate[r0:r1, c0:c1] = out_chunk

        # Temperature Regime
        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                meanT_chunk = self.meanT_daily[r0:r1, c0:c1, :]  # daily temps

                mask_chunk = self.im_mask[r0:r1, c0:c1] != nodata_val if self.set_mask else np.ones((r1 - r0, c1 - c0),
                                                                                                    bool)

                # Compute monthly averages
                meanT_monthly_chunk = np.zeros((r1 - r0, c1 - c0, 12))
                for i in range(r1 - r0):
                    for j in range(c1 - c0):
                        if mask_chunk[i, j]:
                            meanT_monthly_chunk[i, j, :] = averageDailyToMonthly(meanT_chunk[i, j, :], self.leap_year)

                temp_acc_chunk = np.where(meanT_chunk >= 10, meanT_chunk, 0)

                mean_monthly_ge_10 = np.sum(meanT_monthly_chunk >= 10, axis=2)
                mean_monthly_ge_5 = np.sum(meanT_monthly_chunk >= 5, axis=2)
                mean_temp_chunk = np.mean(meanT_chunk, axis=2)
                sum_temp_gt_20 = np.sum(meanT_chunk > 20, axis=2)
                sum_temp_acc = np.sum(temp_acc_chunk, axis=2)

                aez_tclimate_chunk = aez_tclimate[r0:r1, c0:c1]
                temp_conditions = [
                    (mean_monthly_ge_10 == 12) & (mean_temp_chunk >= 20) & mask_chunk,  # TZ1
                    (mean_monthly_ge_5 == 12) & (mean_monthly_ge_10 >= 8) & mask_chunk,  # TZ2
                    (aez_tclimate_chunk == 4) & (mean_monthly_ge_10 >= 5) & (sum_temp_gt_20 >= 75) & (
                                sum_temp_acc > 3000) & mask_chunk,  # TZ3
                    (mean_monthly_ge_10 >= 4) & (mean_temp_chunk >= 0) & mask_chunk,  # TZ4
                    np.isin(mean_monthly_ge_10, [1, 2, 3]) & (mean_temp_chunk >= 0) & mask_chunk,  # TZ5
                    ((mean_monthly_ge_10 == 0) | (mean_temp_chunk < 0)) & mask_chunk  # TZ6
                ]
                temp_values = [1, 2, 3, 4, 5, 6]
                aez_temp_regime[r0:r1, c0:c1] = np.select(temp_conditions, temp_values, default=0)

        # Moisture Regime
        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                lgpt5_chunk = lgpt_5[r0:r1, c0:c1]
                lgp_chunk = lgp[r0:r1, c0:c1]
                lgp_equv_chunk = lgp_equv[r0:r1, c0:c1]

                mask_chunk = self.im_mask[r0:r1, c0:c1] != nodata_val if self.set_mask else np.ones_like(lgpt5_chunk,
                                                                                                         bool)
                out_chunk = np.zeros_like(lgpt5_chunk, dtype=np.int8)

                use_lgp = (lgpt5_chunk > 330) & mask_chunk
                use_lgp_equv = (lgpt5_chunk <= 330) & mask_chunk

                out_chunk = np.where(use_lgp & (lgp_chunk >= 270), 4, out_chunk)
                out_chunk = np.where(use_lgp & (lgp_chunk >= 180) & (lgp_chunk < 270), 3, out_chunk)
                out_chunk = np.where(use_lgp & (lgp_chunk >= 60) & (lgp_chunk < 180), 2, out_chunk)
                out_chunk = np.where(use_lgp & (lgp_chunk >= 0) & (lgp_chunk < 60), 1, out_chunk)

                out_chunk = np.where(use_lgp_equv & (lgp_equv_chunk >= 270), 4, out_chunk)
                out_chunk = np.where(use_lgp_equv & (lgp_equv_chunk >= 180) & (lgp_equv_chunk < 270), 3, out_chunk)
                out_chunk = np.where(use_lgp_equv & (lgp_equv_chunk >= 60) & (lgp_equv_chunk < 180), 2, out_chunk)
                out_chunk = np.where(use_lgp_equv & (lgp_equv_chunk >= 0) & (lgp_equv_chunk < 60), 1, out_chunk)

                aez_moisture_regime[r0:r1, c0:c1] = out_chunk

        # Final AEZ assessment
        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                soil = soil_terrain_lulc[r0:r1, c0:c1]
                temp = aez_temp_regime[r0:r1, c0:c1]
                moist = aez_moisture_regime[r0:r1, c0:c1]

                mask = self.im_mask[r0:r1, c0:c1] != nodata_val if self.set_mask else np.ones_like(soil, bool)

                # Soil based classes
                soil_conditions = [
                    (soil == 8) & mask,
                    (soil == 7) & mask,
                    (soil == 1) & mask,
                    (soil == 6) & mask,
                    (soil == 2) & mask,
                    (soil == 5) & mask
                ]
                soil_values = [56, 57, 49, 51, 52, 50]
                out_soil = np.select(soil_conditions, soil_values, default=0)

                # ---- Desert ----
                desert_mask = (moist == 1) & mask & (out_soil == 0)
                out_desert = np.where(desert_mask, 53, 0)

                # ---- Permafrost ----
                valid_moist = np.isin(moist, [1, 2, 3, 4]) & mask
                out_permafrost = np.select(
                    [
                        (temp == 9) & valid_moist & (out_soil == 0) & (out_desert == 0),
                        (temp == 10) & valid_moist & (out_soil == 0) & (out_desert == 0)
                    ],
                    [54, 55],
                    default=0
                )

                # Formula-based AEZ for remaining cells
                combined_previous = out_soil + out_desert + out_permafrost
                valid_combo = (
                        (temp >= 1) & (temp <= 8) &
                        (moist >= 2) & (moist <= 4) &
                        (soil >= 3) & (soil <= 4) &
                        mask & (combined_previous == 0)
                )
                t = temp - 1
                m = moist - 2
                s = soil - 3
                out_formula = np.zeros_like(out_soil, dtype=np.int8)
                out_formula[valid_combo] = (t[valid_combo] * 6 + m[valid_combo] * 2 + s[valid_combo] + 1)

                # Combining and writing back
                out_chunk = np.where(combined_previous > 0, combined_previous, out_formula)
                aez[r0:r1, c0:c1] = out_chunk


        if self.set_mask:
            return np.where(self.im_mask != nodata_val, aez, np.nan)
        else:
            return aez

    """ 
    Note from Swun: In this code, the logic of temperature amplitude is not added 
    as it brings big discrepency in the temperature regime calculation (in India) 
    compared to previous code. However, the classification schema is now adjusted 
    according to Gunther's agreement and the documentation.
    """
         
    def getMultiCroppingZones(self, t_climate, lgp, lgp_t5, lgp_t10, ts_t10, ts_t0):



        ts_g_t5 = np.zeros((self.im_height, self.im_width), dtype=np.float64)
        ts_g_t10 = np.zeros((self.im_height, self.im_width), dtype=np.float64)

        DAYS_IN_YEAR = self.meanT_daily.shape[2]
        days = np.arange(DAYS_IN_YEAR)
        deg = 5  # polynomial degree

        chunk_size = compute_chunk_size_multi(
            arrays=[ts_g_t5,ts_g_t10],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # interp_daily_temp
                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )

        for r0 in range(0, self.im_height, chunk_size):
            r1 = min(r0 + chunk_size, self.im_height)

            for c0 in range(0, self.im_width, chunk_size):
                c1 = min(c0 + chunk_size, self.im_width)

                # Extract chunk
                temp_chunk = self.meanT_daily[r0:r1, c0:c1, :]  # shape (rows, cols, days)

                # Mask
                if self.set_mask:
                    mask_chunk = self.im_mask[r0:r1, c0:c1] != self.nodata_val
                else:
                    mask_chunk = np.ones((r1 - r0, c1 - c0), dtype=np.bool)

                # Loop over pixels in the chunk
                for i in range(r1 - r0):
                    for j in range(c1 - c0):
                        if not mask_chunk[i, j]:
                            continue

                        temp_1D = temp_chunk[i, j, :]

                        # Polynomial fit
                        polyfit = np.poly1d(np.polyfit(days, temp_1D, deg))
                        interp_daily_temp = polyfit(days)

                        # Vegetative period where temp >= 5
                        veg_period = np.where(interp_daily_temp >= 5)[0]
                        if veg_period.size > 0:
                            start_veg = veg_period[0]
                            end_veg = veg_period[-1] + 1  # +1 to include last day
                        else:
                            start_veg = 0
                            end_veg = DAYS_IN_YEAR

                        # Slice temperature in vegetative period
                        interp_veg = interp_daily_temp[start_veg:end_veg]

                        # Apply thresholds
                        t5 = np.where(interp_veg >= 5, interp_veg, 0)
                        t10 = np.where(interp_veg >= 10, interp_veg, 0)

                        # Sum and store
                        ts_g_t5[r0 + i, c0 + j] = np.sum(t5)
                        ts_g_t10[r0 + i, c0 + j] = np.sum(t10)



        multi_crop_rain = np.zeros((self.im_height, self.im_width), dtype=int)
        nodata_val = self.nodata_val

        for r0 in range(0, self.im_height, chunk_size):
            r1 = min(r0 + chunk_size, self.im_height)

            for c0 in range(0, self.im_width, chunk_size):
                c1 = min(c0 + chunk_size, self.im_width)

                # Extract chunks
                t_climate_chunk = t_climate[r0:r1, c0:c1]
                lgp_chunk = lgp[r0:r1, c0:c1]
                lgp_t5_chunk = lgp_t5[r0:r1, c0:c1]
                lgp_t10_chunk = lgp_t10[r0:r1, c0:c1]
                ts_t0_chunk = ts_t0[r0:r1, c0:c1]
                ts_t10_chunk = ts_t10[r0:r1, c0:c1]
                ts_g_t5_chunk = ts_g_t5[r0:r1, c0:c1]
                ts_g_t10_chunk = ts_g_t10[r0:r1, c0:c1]

                # Mask
                if self.set_mask:
                    mask_chunk = self.im_mask[r0:r1, c0:c1] != nodata_val
                else:
                    mask_chunk = np.ones_like(t_climate_chunk, dtype=np.bool)

                # Initialize chunk output
                out_chunk = np.ones_like(t_climate_chunk, dtype=np.int8)  # default = 1

                # Masks
                mask_tc1 = (t_climate_chunk == 1) & mask_chunk

                cond8 = mask_tc1 & (
                        (lgp_chunk >= 360) & (lgp_t5_chunk >= 360) & (lgp_t10_chunk >= 360) &
                        (ts_t0_chunk >= 7200) & (ts_t10_chunk >= 7000)
                )

                cond6 = mask_tc1 & (
                        (lgp_chunk >= 300) & (lgp_t5_chunk >= 300) & (lgp_t10_chunk >= 240) &
                        (ts_t0_chunk >= 7200) & (ts_g_t5_chunk >= 5100) & (ts_g_t10_chunk >= 4800)
                )

                # Class 4 variants
                cond4a = mask_tc1 & (
                        (lgp_chunk >= 270) & (lgp_t5_chunk >= 270) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 5500) & (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200)
                )

                cond4b = mask_tc1 & (
                        (lgp_chunk >= 240) & (lgp_t5_chunk >= 240) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 6400) & (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200)
                )

                cond4c = mask_tc1 & (
                        (lgp_chunk >= 210) & (lgp_t5_chunk >= 240) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 7200) & (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200)
                )

                cond4 = cond4a | cond4b | cond4c

                # Class 3 variants
                cond3a = mask_tc1 & (
                        (lgp_chunk >= 220) & (lgp_t5_chunk >= 220) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 5500) & (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700)
                )

                cond3b = mask_tc1 & (
                        (lgp_chunk >= 200) & (lgp_t5_chunk >= 200) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 6400) & (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700)
                )

                cond3c = mask_tc1 & (
                        (lgp_chunk >= 180) & (lgp_t5_chunk >= 200) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 7200) & (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700)
                )

                cond3 = cond3a | cond3b | cond3c

                # class 2
                cond2 = mask_tc1 & (
                        (lgp_chunk >= 45) & (lgp_t5_chunk >= 120) & (lgp_t10_chunk >= 90) &
                        (ts_t0_chunk >= 1600) & (ts_t10_chunk >= 1200)
                )

                # Priority Applying
                out_tc1 = np.select(
                    [cond8, cond6, cond4, cond3, cond2],
                    [8, 6, 4, 3, 2],
                    default=1
                )

                mask_tc2 = (t_climate_chunk != 1) & mask_chunk

                cond8_2 = mask_tc2 & (
                        (lgp_chunk >= 360) & (lgp_t5_chunk >= 360) & (lgp_t10_chunk >= 330) &
                        (ts_t0_chunk >= 7200) & (ts_t10_chunk >= 7000)
                )

                cond7_2 = mask_tc2 & (
                        (lgp_chunk >= 330) & (lgp_t5_chunk >= 330) & (lgp_t10_chunk >= 270) &
                        (ts_t0_chunk >= 5700) & (ts_t10_chunk >= 5500)
                )

                cond6_2 = mask_tc2 & (
                        (lgp_chunk >= 300) & (lgp_t5_chunk >= 300) & (lgp_t10_chunk >= 240) &
                        (ts_t0_chunk >= 5400) & (ts_t10_chunk >= 5100) &
                        (ts_g_t5_chunk >= 5100) & (ts_g_t10_chunk >= 4800)
                )

                cond5_2 = mask_tc2 & (
                        (lgp_chunk >= 240) & (lgp_t5_chunk >= 270) & (lgp_t10_chunk >= 180) &
                        (ts_t0_chunk >= 4800) & (ts_t10_chunk >= 4500) &
                        (ts_g_t5_chunk >= 4300) & (ts_g_t10_chunk >= 4000)
                )

                cond4_2 = mask_tc2 & (
                        (lgp_chunk >= 210) & (lgp_t5_chunk >= 240) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 4500) & (ts_t10_chunk >= 3600) &
                        (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200)
                )

                cond3_2 = mask_tc2 & (
                        (lgp_chunk >= 180) & (lgp_t5_chunk >= 200) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 3600) & (ts_t10_chunk >= 3000) &
                        (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700)
                )

                cond2_2 = mask_tc2 & (
                        (lgp_chunk >= 45) & (lgp_t5_chunk >= 120) & (lgp_t10_chunk >= 90) &
                        (ts_t0_chunk >= 1600) & (ts_t10_chunk >= 1200)
                )

                out_tc2 = np.select(
                    [cond8_2, cond7_2, cond6_2, cond5_2, cond4_2, cond3_2, cond2_2],
                    [8, 7, 6, 5, 4, 3, 2],
                    default=1
                )



                out_chunk[mask_tc1] = out_tc1[mask_tc1]
                out_chunk[mask_tc2] = out_tc2[mask_tc2]

                multi_crop_rain[r0:r1, c0:c1] = out_chunk


        multi_crop_irr = np.ones((self.im_height, self.im_width), dtype=np.int8)  # default = 1
        nodata_val = self.nodata_val

        for r0 in range(0, self.im_height, chunk_size):
            r1 = min(r0 + chunk_size, self.im_height)

            for c0 in range(0, self.im_width, chunk_size):
                c1 = min(c0 + chunk_size, self.im_width)

                # Extract chunks
                t_climate_chunk = t_climate[r0:r1, c0:c1]
                lgp_t5_chunk = lgp_t5[r0:r1, c0:c1]
                lgp_t10_chunk = lgp_t10[r0:r1, c0:c1]
                ts_t0_chunk = ts_t0[r0:r1, c0:c1]
                ts_t10_chunk = ts_t10[r0:r1, c0:c1]
                ts_g_t5_chunk = ts_g_t5[r0:r1, c0:c1]
                ts_g_t10_chunk = ts_g_t10[r0:r1, c0:c1]

                # Mask
                if self.set_mask:
                    mask_chunk = self.im_mask[r0:r1, c0:c1] != nodata_val
                else:
                    mask_chunk = np.ones_like(t_climate_chunk, dtype=np.bool)

                # Default chunk output = 1
                out_chunk = np.ones_like(t_climate_chunk, dtype=np.int8)

                # Case 1: t_climate == 1
                mask_tc1 = (t_climate_chunk == 1) & mask_chunk

                out_tc1 = np.select(
                    [
                        mask_tc1 & (lgp_t5_chunk >= 360) & (lgp_t10_chunk >= 360) &
                        (ts_t0_chunk >= 7200) & (ts_t10_chunk >= 7000),

                        mask_tc1 & (lgp_t5_chunk >= 300) & (lgp_t10_chunk >= 240) &
                        (ts_t0_chunk >= 7200) & (ts_g_t5_chunk >= 5100) & (ts_g_t10_chunk >= 4800),

                        mask_tc1 & (lgp_t5_chunk >= 270) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 5500) & (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200),

                        mask_tc1 & (lgp_t5_chunk >= 220) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 5500) & (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700),

                        mask_tc1 & (lgp_t5_chunk >= 120) & (lgp_t10_chunk >= 90) &
                        (ts_t0_chunk >= 1600) & (ts_t10_chunk >= 1200),
                    ],
                    [8, 6, 4, 3, 2],
                    default=1
                ).astype(np.int8)

                #Case 2: t_climate != 1
                mask_tc2 = (t_climate_chunk != 1) & mask_chunk

                out_tc2 = np.select(
                    [
                        mask_tc2 & (lgp_t5_chunk >= 360) & (lgp_t10_chunk >= 330) &
                        (ts_t0_chunk >= 7200) & (ts_t10_chunk >= 7000),

                        mask_tc2 & (lgp_t5_chunk >= 330) & (lgp_t10_chunk >= 270) &
                        (ts_t0_chunk >= 5700) & (ts_t10_chunk >= 5500),

                        mask_tc2 & (lgp_t5_chunk >= 300) & (lgp_t10_chunk >= 240) &
                        (ts_t0_chunk >= 5400) & (ts_t10_chunk >= 5100) &
                        (ts_g_t5_chunk >= 5100) & (ts_g_t10_chunk >= 4800),

                        mask_tc2 & (lgp_t5_chunk >= 270) & (lgp_t10_chunk >= 180) &
                        (ts_t0_chunk >= 4800) & (ts_t10_chunk >= 4500) &
                        (ts_g_t5_chunk >= 4300) & (ts_g_t10_chunk >= 4000),

                        mask_tc2 & (lgp_t5_chunk >= 240) & (lgp_t10_chunk >= 165) &
                        (ts_t0_chunk >= 4500) & (ts_t10_chunk >= 3600) &
                        (ts_g_t5_chunk >= 4000) & (ts_g_t10_chunk >= 3200),

                        mask_tc2 & (lgp_t5_chunk >= 200) & (lgp_t10_chunk >= 120) &
                        (ts_t0_chunk >= 3600) & (ts_t10_chunk >= 3000) &
                        (ts_g_t5_chunk >= 3200) & (ts_g_t10_chunk >= 2700),

                        mask_tc2 & (lgp_t5_chunk >= 120) & (lgp_t10_chunk >= 90) &
                        (ts_t0_chunk >= 1600) & (ts_t10_chunk >= 1200),
                    ],
                    [8, 7, 6, 5, 4, 3, 2],
                    default=1
                ).astype(np.int8)

                out_chunk[mask_tc1] = out_tc1[mask_tc1]
                out_chunk[mask_tc2] = out_tc2[mask_tc2]
                multi_crop_irr[r0:r1, c0:c1] = out_chunk

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
        H, W = self.im_height, self.im_width

        ann_temp_amp = np.zeros((H,W))

        chunk_size = compute_chunk_size_multi(
            arrays=[ann_temp_amp],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # monthly
                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )

        for r0 in range(0,H,chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)
                temp_chunk = self.meanT_daily[r0:r1, c0:c1, :]

                if self.set_mask:
                    mask = self.im_mask[r0:r1, c0:c1] != self.nodata_val
                else:
                    mask = np.ones((r1 - r0, c1 - c0), dtype=np.bool)


                monthly = averageDailyToMonthly(temp_chunk, self.leap_year)

                out_chunk = np.nanmax(monthly, axis=2) - np.nanmin(monthly, axis=2)

                out_chunk[~mask] = 0  # or np.nan
                ann_temp_amp[r0:r1, c0:c1] = out_chunk

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
        if irr_or_rain is None:
            raise Exception('Provide "I" or "R".')

        H, W = self.im_height, self.im_width


        monthly_shortrad = np.zeros((H, W, 12))
        monthly_pr = np.zeros((H, W, 12))

        chunk_size = compute_chunk_size_multi(
            arrays=[monthly_shortrad, monthly_pr],
            extra_arrays=None
        )


        for i0 in range(0, H, chunk_size):
            for j0 in range(0, W, chunk_size):

                i1 = min(i0 + chunk_size, H)
                j1 = min(j0 + chunk_size, W)


                lat_c = self.latitude[i0:i1, j0:j1]
                elev_c = self.elevation[i0:i1, j0:j1]

                minT_c = self.minT_daily[i0:i1, j0:j1, :]
                maxT_c = self.maxT_daily[i0:i1, j0:j1, :]
                shortrad_c = self.shortrad_daily_MJm2day[i0:i1, j0:j1, :]
                wind_c = self.wind_daily[i0:i1, j0:j1, :]
                rh_c = self.rel_humidity_daily[i0:i1, j0:j1, :]
                pr_c = self.totalPrec_daily[i0:i1, j0:j1, :]

                # compute Rn for chunk (parallel inside)
                Rn_chunk = compute_Rn_chunk( #new numba function in ETOCalc that handles parallel processing of chunks
                    lat_c, elev_c,
                    minT_c, maxT_c,
                    shortrad_c, wind_c, rh_c,
                    self.leap_year
                )

                for i in range(i1 - i0):
                    for j in range(j1 - j0):
                        monthly_shortrad[i0 + i, j0 + j, :] = averageDailyToMonthly(
                            Rn_chunk[i, j, :], self.leap_year
                        )
                        monthly_pr[i0 + i, j0 + j, :] = averageDailyToMonthly(
                            pr_c[i, j, :], self.leap_year
                        )


        total_shrad = np.sum(monthly_shortrad, axis=2)
        total_pr = np.sum(monthly_pr, axis=2)

        rdi = np.divide(
            total_shrad,
            total_pr,
            where=total_pr > 0,
            out=np.zeros((H, W))
        )

        if irr_or_rain == 'I':
            rdi = np.minimum(rdi, 1.375)
            npp = np.sum(self.Eto365, axis=2) * rdi * np.exp(-np.sqrt(9.87 + 6.25 * rdi)) * 1000
        else:
            npp = np.sum(self.Eta365, axis=2) * rdi * np.exp(-np.sqrt(9.87 + 6.25 * rdi)) * 1000

        return np.round(npp, 1)


    
    def getBeginningDateofHibernationPeriod(self):
        """
        Calculate the starting date of the hibernation period in a single year frame.
        
        Args:
            None.
        Return:
            lgh_b [2-D NumPy Array]: Beginning date of the hibernation (dormancy period) [Unit: DOY].
        """

        # critical temperature threshold of the least sensitive crop: winter rye
        lgh_b = np.zeros((self.im_height, self.im_width), dtype = np.int8)
        cbtr1, cbtr2 = (-11, -16)

        H, W = self.im_height, self.im_width

        chunk_size = compute_chunk_size_multi(
            arrays=[lgh_b],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # monthly
                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )

        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                chunk = self.meanT_daily[r0:r1, c0:c1, :]
                if self.set_mask:
                    mask = self.im_mask[r0:r1, c0:c1] != self.nodata_val
                else:
                    mask = np.ones((r1 - r0, c1 - c0), dtype=np.bool)


                monthly = averageDailyToMonthly(chunk, leap_year)
                tadif0 = np.nanmax(monthly, axis=2) - np.nanmin(monthly, axis=2)

                cond1 = (tadif0 < 35) & (tadif0 > 20) & mask
                cond2 = (tadif0 > 35) & mask
                cond3 = (tadif0 < 20) & mask

                out_chunk = np.zeros((r1 - r0, c1 - c0), dtype=np.float32)


                out_chunk[cond1] = cbtr1
                out_chunk[cond2] = cbtr1 + (cbtr2 - cbtr1) * (35 - tadif0) / 15
                out_chunk[cond3] = cbtr2

                cbtr = out_chunk[:, :, np.newaxis]

                # three criteria must be satisfied for dormancy (hibernation period) determination
                # 1. Average Temperature must be less than 5
                # 2. Dormancy period must be less than 200 days
                # 3. Average Temperature must satisfy crop-specific temperature threshold.

                dormancy_days = (chunk < 5) & (chunk >= cbtr)

                idx = np.argmax(dormancy_days, axis=2)

                # handle no True case
                no_dormancy = ~np.any(dormancy_days, axis=2)

                lgh_b_chunk = idx + 1
                lgh_b_chunk[no_dormancy] = 1

                # assign back
                lgh_b[r0:r1, c0:c1] = lgh_b_chunk

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

        lgh = np.zeros((self.im_height, self.im_width), dtype=np.int8)
        cbtr1, cbtr2 = (-11, -16)

        H, W = self.im_height, self.im_width

        chunk_size = compute_chunk_size_multi(
            arrays=[lgh],
            extra_arrays=[
                # Monthly arrays (dominant)
                {"shape_factor": (12,), "dtype": np.float64},  # monthly
                # Output
                {"shape_factor": (1,), "dtype": np.int8},

                # Masks (approximate)
                {"shape_factor": (1,), "dtype": np.bool_},
            ],
        )

        for r0 in range(0, H, chunk_size):
            r1 = min(r0 + chunk_size, H)
            for c0 in range(0, W, chunk_size):
                c1 = min(c0 + chunk_size, W)

                chunk = self.meanT_daily[r0:r1, c0:c1, :]
                if self.set_mask:
                    mask = self.im_mask[r0:r1, c0:c1] != self.nodata_val
                else:
                    mask = np.ones((r1 - r0, c1 - c0), dtype=np.bool)

                monthly = averageDailyToMonthly(chunk, leap_year)
                tadif0 = np.nanmax(monthly, axis=2) - np.nanmin(monthly, axis=2)

                cond1 = (tadif0 < 35) & (tadif0 > 20) & mask
                cond2 = (tadif0 > 35) & mask
                cond3 = (tadif0 < 20) & mask

                out_chunk[cond1] = cbtr1
                out_chunk[cond2] = cbtr1 + (cbtr2 - cbtr1) * (35 - tadif0) / 15
                out_chunk[cond3] = cbtr2

                dormancy_days = (temp_chunk < 5) & (temp_chunk >= cbtr)

                # apply mask
                dormancy_days &= mask[:, :, None]

                # duration (count of True days)
                lgh_chunk = np.sum(dormancy_days, axis=2)

                # assign back
                lgh[r0:r1, c0:c1] = lgh_chunk
        return lgh


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


