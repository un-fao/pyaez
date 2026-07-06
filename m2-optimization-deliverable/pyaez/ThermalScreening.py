"""
PyAEZ version 2.4 (Dec 2024)
Thermal Screening
2020: N. Lakmal Deshapriya
2022/2023: Swun Wunna Htet
2023 (Dec): Swun Wunna Htet
2024 (Dec): Swun Wunna Htet

Modification:
1. LGPT screening is now removed because of the new consideration with LGPT and LGP in Module 2.
2. Smoothening is applied to TSUM Screening based on GAEZ routine. However, the investigation with
   residual handling must be asserted.
3. Object classes will be used for input data provision. The calculation 
   procedures will be omitted out from object class for possible Numba enhancement (but not every single functions).
"""

import numpy as np
import numba as nb
from numba.typed import List
from scipy.interpolate import interp1d, PchipInterpolator  # shape-preserving cubic
from typing import NamedTuple

np.round_ = np.round


class CropSpecificYearContext(NamedTuple):
    """Year-constant inputs for crop-specific temperature profile screening."""
    rule: object
    constr_type: list
    optimal: list
    sub_optimal: list
    not_suitable: list
    tpro365: tuple
    RHavg: float
    RHmin: float
    RHmax: float
    LGPT00: object
    LGPT05: object
    LGPT10: object

#----------------------------------------------Major Functions Starts Here----------------------------------------------------------------
# This numba enhancement decorator is omitted out because the function runs faster without nb.jit() (nb.jit takes 2.3 seconds, no decorator produce faster 0.5 seconds)
# @nb.jit(nopython = True)
# perennial flag in function removed as it is not applied.
def getReductionFactorNumba(set_Tsum_screening:bool, LnS, LsO, LO, HnS, HsO, HO, tsum0,
                            set_CropSpecificRule:bool, crop_specific_rule_data, hibernating_flag, vern_factor):
    """Calculation of LUT specific thermal suitability factor (fc1)
    
    Args:
        set_Tsum_screening (Bool): TSUM screening activation
        LnS (int): Lower boundary of not-suitable accumulated heat unit range.
        LsO (int): Lower boundary of sub-optimal accumulated heat unit range.
        LO (int): Lower boundary of optimal accumulated heat unit range.
        HnS (int):Upper boundary of not-suitable accumulated heat range.
        HsO (int): Upper boundary of sub-optimal accumulated heat range.
        HO (int): Upper boundary of not-suitable accumulated heat range.
        tsum0 (float): Temperature summation at 0 Deg Threshold
        set_CropSpecificRule (bool): Crop-specific Rule activation
        crop_specific_rule_data (list): calculated crop-specific rule constraints
        hibernating_flag (bool): Hibernation activation
        vern_factor (float): calculated vernalization factor.
        
    Return:
        fc1 (float): therml suitability factor, value between 0 (not suitable) and 1  (suitable)
    """
    
    fc1_final = 1.

    # TSUM screening
    if set_Tsum_screening:
        # Start TSUM screening
        # Optimal range
        if tsum0 in range(LO, HO):
            f1 = 1.
            fc1_final = min(f1, fc1_final)

        # Within Sub-optimal range (Part 1) (25% reduction factor)
        elif tsum0 in range(LsO, LO):
            f1 = ((tsum0-LsO)/(LO-LsO)) * 0.25 + 0.75
            fc1_final = min(f1, fc1_final)

        # Within Sub-optimal range (Part 2) (25% reduction factor)
        elif tsum0 in range(HO, HsO):
            f1 = ((HsO-tsum0)/(HsO-HO)) * 0.25 + 0.75
            fc1_final = min(f1, fc1_final)

        # Within Marginal range (Part 1) (75% reduction factor)
        elif tsum0 in range(LnS, LsO):
            f1 = ((tsum0-LnS)/(LsO-LnS)) * 0.75
            fc1_final = min(f1, fc1_final)

        # Within Marginal range (Part 2) (75% reduction factor)
        elif tsum0 in range(HsO, HnS):
            f1 = ((HnS-tsum0)/(HnS-HsO)) * 0.75
            fc1_final = min(f1, fc1_final)

        # Within Not suitable range (100% reduction factor)
        elif tsum0 <= LnS or tsum0 >= HnS:
            f1 = 0
            fc1_final = min(f1, fc1_final)

    # Vernalization Screening
    if hibernating_flag:
        fc1_final = min(vern_factor, fc1_final)

    # Crop Specific Rule Screening (Temperature Profile Screening)
    if set_CropSpecificRule:
        
        specific_data = crop_specific_rule_data
        calc_value = specific_data[0]
        constr_type =  specific_data[1]
        optimal =  specific_data[2]
        sub_optimal =  specific_data[3]
        not_suitable =  specific_data[4]

        for i in range(len(calc_value)):
            # """Check if all threshold values are the same"""
            if constr_type[i] == '<=' or constr_type[i] == '≤':
                # """Calculated value will be compared with optimum threshold"""
                if optimal[i] == sub_optimal[i] == not_suitable[i]:
                    if calc_value[i] <= optimal[i]:
                        f1 = vern_factor if hibernating_flag else 1
                    else:
                        f1 = 0
                    fc1_final = min(f1, fc1_final)
                
                elif optimal[i] != sub_optimal[i] == not_suitable[i]:
                    if calc_value[i] <= optimal[i]:
                        f1 = vern_factor if hibernating_flag else 1
                    # """If calculated value within range between optimal and sub-optimum/not-suitable"""
                    elif calc_value[i] > optimal[i] and calc_value[i] <= sub_optimal[i]:
                        f1 = vern_factor if hibernating_flag else ((calc_value[i] - optimal[i])/(sub_optimal[i] - optimal[i]) * 0.25) + 0.75
                    else:
                        f1 = 0
                    fc1_final = min(f1, fc1_final)

                 # """If all thresholds are different, go linear interpolation to each threshold interval"""
                elif optimal[i] != sub_optimal[i] != not_suitable[i]:
                    if calc_value[i] <= optimal[i]:
                        f1 = vern_factor if hibernating_flag else 1
                    # """If calculated value within range between optimal and sub-optimum/not-suitable"""
                    elif calc_value[i] > optimal[i] and calc_value[i] <= sub_optimal[i]:
                        f1 = ((calc_value[i] - optimal[i])/(sub_optimal[i] - optimal[i]) * 0.25) + 0.75
                    # """For calculated value beyond sub-optimum/not-suitable, use previous linear interpolation (But not sure)"""
                    elif calc_value[i] > sub_optimal[i] and calc_value[i] <= not_suitable[i]:
                        f1 = ((calc_value[i] - not_suitable[i])/(sub_optimal[i] - not_suitable[i]) * 0.25) + 0.75
                    else:
                        f1 = 0
                    fc1_final = min(f1, fc1_final)
            
            # """Check if all threshold values are the same"""
            elif constr_type[i] == '>=' or constr_type[i] == '≥':
                if optimal[i] == sub_optimal[i] == not_suitable[i]:
                    f1 = 1 if calc_value[i] >= optimal[i] else 0
                    fc1_final = min(f1, fc1_final)
                elif optimal[i] != sub_optimal[i] == not_suitable[i]:
                    if calc_value[i] >= optimal[i]:
                        f1 = vern_factor if hibernating_flag else 1
                    elif calc_value[i] < optimal[i] and calc_value[i] >= sub_optimal[i]:
                        f1 = vern_factor if hibernating_flag else ((calc_value[i] - optimal[i])/(sub_optimal[i] - optimal[i]) * 0.25) + 0.75
                    else:
                        f1 = 0
                    fc1_final = min(f1, fc1_final)
                
                elif optimal[i] != sub_optimal[i] != not_suitable[i]:
                    if calc_value[i] >= optimal[i]:
                        f1 = vern_factor if hibernating_flag else 1
                    elif calc_value[i] < optimal[i] and calc_value[i] >= sub_optimal[i]:
                        f1 = ((calc_value[i] - optimal[i])/(sub_optimal[i] - optimal[i]) * 0.25) + 0.75
                    elif calc_value[i] < sub_optimal[i] and calc_value[i] >= not_suitable[i]:
                        f1 = ((calc_value[i] - not_suitable[i])/(sub_optimal[i] - not_suitable[i]) * 0.25) + 0.75
                    else:
                        f1 = 0
                    fc1_final = min(f1, fc1_final)

    return fc1_final


# ---------------------------------------- Numba Enhanced Functions End Here---------------------------------- #
# ------------------Intermediate Functions (Not available for Numba enhancement) Starts Here -------------------#
def getInterpolatedTempData(temp1D, interp_method = 'Pchip'):
    """
    Get smoothened temperature data for 365 days.
    
    Parameter
    ---------
    temp1D (1-D NumPy Array): Input mean temperature (Deg C)
    interp_method (string): interpolation methodology, can be either 'interp1' or 'Pchip' (default)
    
    Returns
    -------
    interp_temp1D (1-D NumPy Array): Smoothened mean temperature.
    """

    if interp_method == 'interp1':
        interp1D = np.zeros(temp1D.shape)
    
        tmp = np.array([temp1D[k-1] for k in range(15,temp1D.shape[0]+1,30)])
        mid_doy = np.arange(15,temp1D.shape[0]+1,30)
    
        # Quadratic spline interpolation for 330 days
        int_mdl = interp1d(mid_doy, tmp, kind='quadratic', fill_value='extrapolate')
    
        interp1D = int_mdl(np.arange(1,temp1D.shape[0]+1))

    elif interp_method == 'Pchip':
        n = temp1D.shape[0]
        mid_doy = np.arange(15, n+1, 30)
        anchors_x = np.unique(np.r_[1, mid_doy, n])
        tmp = np.array([temp1D[k-1] for k in anchors_x])  # same as temp1D[anchors_x-1]
        
        int_mdl = PchipInterpolator(anchors_x, tmp, extrapolate=True)
        interp1D = int_mdl(np.arange(1, n+1))
        
        # Enforce bounds: choose either the empirical [min,max] or a known domain
        lo, hi = float(np.nanmin(tmp)), float(np.nanmax(tmp))
        # If you know temps must be, say, [0, 1], use: lo, hi = 0.0, 1.0
        interp1D = np.clip(interp1D, lo, hi)

    return interp1D

def findbeginningLGPT(mean_temp, threshold):
    """
    Find the beginning date of LGPt with specific threshold
    
    Parameter
    ---------
    mean_temp (1-D NumPy Array): smoothed mean temperature [Deg C]
    threshold (float/int)

    Return
    ------
    lgbt (int): beginning date of the threshold threshold
    """
    for i in range(mean_temp.shape[0]):
        if mean_temp[i]>threshold:
            return i
        else:
            continue
    return i

def getTemperatureGrowingPeriod(temp1D, threshold):
    """
    Calculation of thermal growing period at defined threshold.
    
    Parameters
    ----------
    temp1D (1-D NumPy Array): Input mean temperature (Deg C)
    
    Returns
    -------
    LGPt (int): thermal growing period.
    """
    # Calculation of Temp Profile for 1-D numpy array of climate data input
    interp1D = getInterpolatedTempData(temp1D)

    lgpt = interp1D>=threshold

    return np.nansum(lgpt)

def getTemperatureSum(temp1D, threshold):
    """
    Calculation of temperature summation at defined threshold.
    
    Parameters
    ----------
    temp1D (1-D NumPy Array): Input mean temperature (Deg C)
    
    Returns
    -------
    TSUM (int): temperature summation.
    """
    # Calculation of Temp Profile for 1-D numpy array of climate data input
    interp1D = getInterpolatedTempData(temp1D)

    interp1D[interp1D <= threshold] = 0

    return np.round(np.sum(interp1D), decimals=0)

def getSmoothTemp(temp1D):

    """Get a smoothened temperature curve done by quadratic spline.
    For 365/366 days.
    
    Parameters
    ----------
    temp1D (1-D NumPy Array): Input mean temperature (Deg C)
    
    Return
    ------
    smootheTemp (1-D NumPy Array): Quadratic spline smoothened temperature (Deg C)
    """
    # Calculation of Temp Profile for 1-D numpy array of climate data input
    interp1D = getInterpolatedTempData(temp1D)

    return interp1D


def getTempTrend(temp1D):
    """
    Create a smoothened temperature curve of 365/366 days.
    
    Parameters
    ----------
    temp1 (1-D NumPy Array): Input mean temperature (Deg C)
    
    Returns
    -------
    Temperature Trend (1-D NumPy Array): Upward (+1)/ Downward (-1) trend data.
    """
    # Calculation of Temp Profile for 1-D numpy array of climate data input
    interp1D = getInterpolatedTempData(temp1D)
    
    # Detect the warmest and coldest day of year.
    tmaxidx = np.argmax(interp1D)
    tminidx= np.argmin(interp1D)

    # Any days between the interval of warmest and coldest DOY will be decreasing trend
    meanT_diff = np.ones(interp1D.shape)
    meanT_diff[tmaxidx:tminidx] = -1

    return meanT_diff

def getTemperatureProfile(temp1D):
    """
    Calculation of temperature profile. The length of temperature data differs depend on 
    crop type (annuals or perennials).
    
    Parameters
    ----------
    temp1 (1-D NumPy Array): Input mean temperature (Deg C)
    
    Returns
    -------
    None.
    """
    # Calculation of Temp Profile for 1-D numpy array of climate data input
    interp1D = getInterpolatedTempData(temp1D)
    
    # Detect the warmest and coldest day of year.
    tmaxidx = np.argmax(interp1D)
    tminidx= np.argmin(interp1D)

    # Any days between the interval of warmest and coldest DOY will be decreasing trend
    meanT_diff = np.ones(interp1D.shape)
    meanT_diff[tmaxidx:tminidx] = -1

    # Allocating the temperature to the corresponding classes. New clases introduced (A0, B0)
    A9 = np.sum(np.logical_and(meanT_diff > 0, interp1D < -5))
    A8 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= -5, interp1D < 0)))
    A7 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 0, interp1D < 5)))
    A6 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 5, interp1D < 10)))
    A5 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 10, interp1D < 15)))
    A4 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 15, interp1D < 20)))
    A3 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 20, interp1D < 25)))
    A2 = np.sum(np.logical_and(meanT_diff > 0, np.logical_and(
        interp1D >= 25, interp1D < 30)))
    A1 = np.sum(np.logical_and(meanT_diff > 0, interp1D >= 30, interp1D < 35))

    A0 = np.sum(np.logical_and(meanT_diff >0, interp1D >=35))

    B9 = np.sum(np.logical_and(meanT_diff < 0, interp1D < -5))
    B8 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= -5, interp1D < 0)))
    B7 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 0, interp1D < 5)))
    B6 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 5, interp1D < 10)))
    B5 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 10, interp1D < 15)))
    B4 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 15, interp1D < 20)))
    B3 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 20, interp1D < 25)))
    B2 = np.sum(np.logical_and(meanT_diff < 0, np.logical_and(
        interp1D >= 25, interp1D < 30)))
    B1 = np.sum(np.logical_and(meanT_diff < 0, interp1D >= 30, interp1D < 35))
    B0 = np.sum(np.logical_and(meanT_diff <0, interp1D >=35))

    # releasing memory
    del (temp1D, interp1D, meanT_diff)

    return [A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9]


def insertCropSpecificRuleParameters(data):

    rule = data['Constraint'].to_numpy()
    constr_type = data['Type'].to_numpy()
    optimal = list(data['Optimal'].to_numpy())
    sub_optimal = list(data['Sub-Optimal'].to_numpy())
    not_suitable = list(data['Not-Suitable'].to_numpy())

    return rule, constr_type, optimal, sub_optimal, not_suitable


def _eval_crop_specific_rule_expressions(
        rule,
        tpro365,
        trpocycle,
        RHavg: float,
        RHmin: float,
        RHmax: float,
        LGPT00,
        LGPT05,
        LGPT10,
    ):
    """Evaluate crop-specific constraint strings (same local namespace as legacy path)."""
    N0a = tpro365[0]
    N1a = tpro365[1]
    N2a = tpro365[2]
    N3a = tpro365[3]
    N4a = tpro365[4]
    N5a = tpro365[5]
    N6a = tpro365[6]
    N7a = tpro365[7]
    N8a = tpro365[8]
    N9a = tpro365[9]
    N0b = tpro365[10]
    N1b = tpro365[11]
    N2b = tpro365[12]
    N3b = tpro365[13]
    N4b = tpro365[14]
    N5b = tpro365[15]
    N6b = tpro365[16]
    N7b = tpro365[17]
    N8b = tpro365[18]
    N9b = tpro365[19]
    N0 = N0a + N0b
    N1 = N1a + N1b
    N2 = N2a + N2b
    N3 = N3a + N3b
    N4 = N4a + N4b
    N5 = N5a + N5b
    N6 = N6a + N6b
    N7 = N7a + N7b
    N8 = N8a + N8b
    N9 = N9a + N9b

    L0a = trpocycle[0]
    L1a = trpocycle[1]
    L2a = trpocycle[2]
    L3a = trpocycle[3]
    L4a = trpocycle[4]
    L5a = trpocycle[5]
    L6a = trpocycle[6]
    L7a = trpocycle[7]
    L8a = trpocycle[8]
    L9a = trpocycle[9]
    L0b = trpocycle[10]
    L1b = trpocycle[11]
    L2b = trpocycle[12]
    L3b = trpocycle[13]
    L4b = trpocycle[14]
    L5b = trpocycle[15]
    L6b = trpocycle[16]
    L7b = trpocycle[17]
    L8b = trpocycle[18]
    L9b = trpocycle[19]
    L0 = L0a + L0b
    L1 = L1a + L1b
    L2 = L2a + L2b
    L3 = L3a + L3b
    L4 = L4a + L4b
    L5 = L5a + L5b
    L6 = L6a + L6b
    L7 = L7a + L7b
    L8 = L8a + L8b
    L9 = L9a + L9b

    calc_value = []
    for i in range(len(rule)):
        calc_value.append(eval(rule[i]))
    return calc_value


def build_crop_specific_year_context(
        data,
        input_temp,
        input_RH,
        LGPT00,
        LGPT05,
        LGPT10,
    ) -> CropSpecificYearContext:
    """
    Precompute year-level temperature profile bins and rule metadata for one pixel.
    Safe to reuse across all planting-date cycles at that location.
    """
    rule, constr_type, optimal, sub_optimal, not_suitable = insertCropSpecificRuleParameters(data)
    tpro365 = tuple(getTemperatureProfile(input_temp))
    # Legacy unpack order matches calculateTemperatureProfileClasses (array is min, max, mean).
    _rh = getRHstats(input_RH)
    RHavg, RHmin, RHmax = _rh[0], _rh[1], _rh[2]
    return CropSpecificYearContext(
        rule=rule,
        constr_type=constr_type,
        optimal=optimal,
        sub_optimal=sub_optimal,
        not_suitable=not_suitable,
        tpro365=tpro365,
        RHavg=RHavg,
        RHmin=RHmin,
        RHmax=RHmax,
        LGPT00=LGPT00,
        LGPT05=LGPT05,
        LGPT10=LGPT10,
    )


def calculateTemperatureProfileClasses_for_cycle(
        year_ctx: CropSpecificYearContext,
        input_temp,
        i_cycle: int,
        cycle_len: int,
    ):
    """Cycle-window portion of calculateTemperatureProfileClasses (uses cached year context)."""
    trpocycle = getTemperatureProfile(input_temp[i_cycle:i_cycle + cycle_len])
    calc_value = _eval_crop_specific_rule_expressions(
        year_ctx.rule,
        year_ctx.tpro365,
        trpocycle,
        year_ctx.RHavg,
        year_ctx.RHmin,
        year_ctx.RHmax,
        year_ctx.LGPT00,
        year_ctx.LGPT05,
        year_ctx.LGPT10,
    )
    return (
        calc_value,
        year_ctx.constr_type,
        year_ctx.optimal,
        year_ctx.sub_optimal,
        year_ctx.not_suitable,
    )


def getRHstats(RelHum_1d: np.ndarray, ax=None, assume_percent=False, debug=False):
    """
    Compute stats and plot RH using a monotone, shape-preserving interpolation (PCHIP).
    Values are clipped to [0, 1].
    """
    RelHum_1d = np.asarray(RelHum_1d, dtype=float)
    n = RelHum_1d.shape[0]
    if n == 0:
        raise ValueError("RelHum_1d is empty.")

    # If the data are in percent [0..100], convert to fraction [0..1]
    if assume_percent or np.nanmax(RelHum_1d) > 1.5:
        RelHum_1d = RelHum_1d / 100.0

    # Choose anchors: use mid-of-month, but also include endpoints to stabilize edges
    mid_doy = np.arange(15, n + 1, 30)           # 1-based days: 15, 45, 75, ...
    anchors_x = np.unique(np.r_[1, mid_doy, n])  # add day 1 and day n, ensure unique
    anchors_y = RelHum_1d[anchors_x - 1]         # zero-based indexing

    # Build monotone cubic interpolator (no crazy overshoot)
    mdl = PchipInterpolator(anchors_x, anchors_y, extrapolate=True)

    days = np.arange(1, n + 1)
    interp1D = mdl(days)

    # Clip to [0, 1] because RH (fraction) must be within that interval
    interp1D = np.clip(interp1D, 0.0, 1.0)

    # Plot
    if debug:
        if ax is None:
            ax = plt.gca()
        ax.plot(days, interp1D, color='tab:blue', lw=2, label='Interpolated RH (PCHIP)')
        ax.scatter(days, RelHum_1d, s=22, color='tab:orange', alpha=0.85, label='Original RH (daily)')
        ax.scatter(anchors_x, anchors_y, s=40, color='tab:green', edgecolor='black', zorder=3,
                   label='Anchor points')
        ax.set_xlabel('Day')
        ax.set_ylabel('Relative Humidity (fraction)')
        ax.set_title('Daily RH: Original Points and Shape-Preserving Interpolation')
        ax.grid(True, alpha=0.3)
        ax.legend()

    RHmin = float(np.nanmin(interp1D))
    RHmax = float(np.nanmax(interp1D))
    RHavg = float(np.nanmean(interp1D))
    return np.array([RHmin, RHmax, RHavg], dtype=float)
    
# 4 Modification

def calculateTemperatureProfileClasses(
        data,
        input_temp,
        input_RH,
        LGPT00,
        LGPT05,
        LGPT10,
        i_cycle,
        cycle_len
    ):
    """
    Compute temperature‑profile class metrics for a full year (365 days) and for a
    crop cycle window, then evaluate user/crop‑specific constraint equations.

    The function:
      1) Loads crop‑specific rule parameters via `insertCropSpecificRuleParameters(data)`,
      2) Computes temperature‑profile counts for the whole year (`tpro365`) and for the
         selected cycle window (`trpocycle`) using `getTemperatureProfile(...)`,
      3) Aggregates paired daytime/nighttime (or "a"/"b") bins into 10 combined bins
         for the year (N0..N9) and for the cycle (L0..L9),
      4) Computes cycle‑window RH statistics (`RHavg`, `RHmin`, `RHmax`) from `input_RH`,
      5) Evaluates each constraint expression in `rule` (strings) using the variables
         defined above, returning the numeric values and the associated class thresholds.

    Parameters
    ----------
    data : Any
        Crop‑specific configuration or rule source passed to
        `insertCropSpecificRuleParameters(data)`. The helper is expected to return a
        5‑tuple/list: `(rule, constr_type, optimal, sub_optimal, not_suitable)`.
        - `rule` is a sequence of string expressions to be evaluated (e.g., "N3 >= 10 and L4 > 5").
        - `constr_type`, `optimal`, `sub_optimal`, `not_suitable` are parallel sequences
          describing constraint labels and class thresholds.
    input_temp : array_like
        1‑D sequence of daily mean (or representative) temperatures for the full year.
        Expected length: 365 (indexable). Units should match what
        `getTemperatureProfile(...)` expects (e.g., °C).
    input_RH : array_like
        1‑D sequence of daily relative humidity values (percent, 0–100) for the full year.
        Used to compute `RHavg`, `RHmin`, `RHmax` over the crop cycle window.
    LGPT00 : int
        LGPt0 value for the specific crop
    LGPT05 : int
        LGPt5 value for the specific crop
    LGPT10 : int
        LGPt10 value for the specific crop
    i_cycle : int
        Start index (0‑based) of the crop cycle window within the year. The slice
        `[i_cycle : i_cycle + cycle_len]` must lie within the bounds of `input_temp`
        and `input_RH`.
    cycle_len : int
        Length (in days) of the crop cycle window to analyze.

    Returns
    -------
    calc_value : list of numbers
        Numeric results of evaluating each rule expression in `rule`. Each element is
        typically a boolean (0/1) or numeric score, depending on the expressions.
    constr_type : list
        Constraint identifiers (e.g., names or codes) corresponding to each rule.
    optimal : list
        Thresholds/parameters defining **optimal** class for each constraint.
    sub_optimal : list
        Thresholds/parameters defining **sub‑optimal** class for each constraint.
    not_suitable : list
        Thresholds/parameters defining **not‑suitable** class for each constraint.

    Notes
    -----
    - **Temperature profiles**: `getTemperatureProfile(x)` is expected to return a 20‑element
      vector where indices 0–9 (`*a`) and 10–19 (`*b`) are paired bins; the function
      aggregates them as:
        N*k = N*k_a + N*k_b  (yearwide)  and  L*k = L*k_a + L*k_b (cycle window),
      yielding k ∈ {0..9}. The exact bin meaning (e.g., temperature ranges) depends on
      your `getTemperatureProfile` implementation.
    - **Relative humidity**: `getRHstats(...)` must return `(RHavg, RHmin, RHmax)` for the
      provided window; units should match the rules’ expectations.
    - **Variables available to rules**:
        * Year bins: `N0..N9`
        * Cycle bins: `L0..L9`
        * Cycle RH stats: `RHavg`, `RHmin`, `RHmax`
      If your rules reference other variables (e.g., `LGPT00`), ensure they are defined
      before rule evaluation.
    - **Safety**: Rule expressions are evaluated with Python `eval`. Only use trusted rule
      sources or replace `eval` with a safe expression parser if rules can be user‑provided.
    - **Bounds**: The slice `[i_cycle : i_cycle + cycle_len]` must be valid for both
      `input_temp` and `input_RH`. Validate indices before calling if they may vary.

    Examples
    --------
    >>> # Minimal illustration (assuming helpers are available in scope):
    >>> calc_value, ctype, opt, subopt, ns = calculateTemperatureProfileClasses(
    ...     data=cfg,                       # source for insertCropSpecificRuleParameters
    ...     input_temp=temp365,             # length 365
    ...     input_RH=rh365,                 # length 365
    ...     LGPT00=None,                    # lgpt0
    ...     LGPT05=None,                    # lgpt5
    ...     LGPT10=None,                    # lgpt10
    ...     i_cycle=90,                     # start at day 90
    ...     cycle_len=120                   # 120-day crop    ...     cycle_len=120                   # 120-day crop cycle
    ... )
    """
    
    year_ctx = build_crop_specific_year_context(
        data, input_temp, input_RH, LGPT00, LGPT05, LGPT10
    )
    return calculateTemperatureProfileClasses_for_cycle(
        year_ctx, input_temp, i_cycle, cycle_len
    )

#----------------------------------------------Major Functions Ends Here----------------------------------------------------------------
