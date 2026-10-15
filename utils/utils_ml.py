import xdas.signal as xp
import pandas as pd
import numpy as np
import xarray as xr

from pyproj import Geod

def to_timedelta64(t):
    return np.rint(t * 1e9).astype("timedelta64[ns]")


def nansmad(x):
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def window_timing(picks_df, da, pre_win=0.2, wlen=8):
    """
    Funtion to generate datasets of winow start and end times and actual picks.
    Input:
        picks_df: pandas dataframe with picks
        da: xdas dadatarray with distance and time
        pre_win: pre-window time in seconds
        wlen: window length in seconds
    Output:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        toa: xarray dataset with actual picks for each phase
    """
    toa = xr.Dataset()
    starttime = xr.Dataset()
    endtime = xr.Dataset()
    for phase in picks_df["phase"].unique():
        phase_df = picks_df[picks_df["phase"] == phase]
        phase_df = phase_df.drop_duplicates(subset=["offset"])
        phase_df = phase_df.sort_values("offset")
        da_rounded = np.round(da["distance"].values, 3)
        phase_df_rounded = phase_df.copy()
        phase_df_rounded["offset"] = np.round(phase_df_rounded["offset"].values, 3)
        da_df = pd.DataFrame({"distance": da_rounded})
        merged = da_df.merge(phase_df_rounded[["offset", "time"]], left_on="distance", right_on="offset", how="left")
        toa[phase] = xr.DataArray(merged["time"].values, {"distance": da["distance"].values})
        starttime[phase] = xr.DataArray(merged["time"].values - to_timedelta64(pre_win), {"distance": da["distance"].values})
        endtime[phase] = xr.DataArray(merged["time"].values + to_timedelta64(wlen - pre_win), {"distance": da["distance"].values})
    return starttime, endtime, toa


def clean_picks_smad(window_timing, threshold=10):
    """
    Function to clean picks using a smad filter.
    Input:
        window_timing: dataset with window times for each phase
        threshold: threshold for smad filter
    Output:
        cleaned dataset
    """
    for phase in window_timing:
        picks = (window_timing[phase].values - np.datetime64(0, "s")) / np.timedelta64(1, "s")
        picks = picks - np.nanmedian(picks)
        nsmad = (picks - nansmad(picks)) / nansmad(picks)
        mask = np.abs(nsmad) <= threshold
        window_timing[phase].values[~mask] = np.datetime64("NaT")
    return window_timing


def add_noise(starttime, endtime, da, noise_wlen):
    """
    Function to add noise start and end times.
    Input:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        da: xdas dadatarray with distance and time
        noise_wlen: noise time length time in seconds
    Output:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
    """
    starttime["Noise"] = starttime["P"].copy(
        data=np.full(
            len(starttime["P"]),
            np.min(da["time"].values),
            dtype="datetime64[ns]",
        )
    )
    endtime["Noise"] = endtime["P"].copy(
        data=starttime["Noise"] + to_timedelta64(noise_wlen)
    )
    return starttime, endtime


def avoid_overlap(starttime, endtime, phase1, phase2):
    """
    Function to avoid overlap between two phases.
    Input:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        phase1: first phase
        phase2: second phase
    Output: 
        starttime: cleaned dataset after overlapping removal
        endtime: cleaned dataset after overlapping removal
    """
    safe_delay = (starttime[phase2] - starttime[phase1]) / 8
    mask = endtime[phase1] >= starttime[phase2] - safe_delay
    endtime[phase1][mask] = (
        starttime[phase2][mask] - safe_delay[mask]
    )
    return starttime, endtime


def correct_pwin(starttime, endtime, wlen):
    """
    Function to correct window length when successive phases do  not have picks
    Input:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        wlen: window length in seconds
    Output:
        starttime: cleaned dataset after window length correction
        endtime: cleaned dataset after window length correction
    """ 
    mask = np.isnat(starttime["S"])
    wlen_p = np.nanpercentile(starttime["S"].values - starttime["P"].values, 10)
    if wlen_p > to_timedelta64(wlen):
        wlen_p = to_timedelta64(wlen)
    endtime["P"].values[mask] = starttime["P"].values[mask] + wlen_p
    return starttime, endtime


def avoid_nullwin(starttime, endtime, phase):
    """
    Function to avoid null windows.
    Input:
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        phase: phase to check for null windows
    Output:
        starttime: cleaned dataset after null window removal
        endtime: cleaned dataset after null window removal
    """
    mask = starttime[phase] >= endtime[phase]
    starttime[phase].values[mask] = np.datetime64("NaT", "ns")
    endtime[phase].values[mask] = np.datetime64("NaT", "ns")
    return starttime, endtime


def preprocess_disp(da, tapering_frac, lowcut, highcut, avg_len):
    """
    Function to preprocess displacement data.
    Input:
        da: xdas dadatarray with distance and time
        tapering_frac: fraction of tapering window
        lowcut: low cut frequency in Hz
        highcut: high cut frequency in Hz
        avg_len: window length for m moving average in meters
    Output:
        da: preprocessed displacement dataarray
    """
    da = da - da.mean(dim="distance")
    da = xp.integrate(da, dim="distance")
    da = xp.sliding_mean_removal(da, wlen=avg_len, dim="distance")
    da = xp.integrate(da, dim="time")
    da = xp.taper(da, window=("tukey", tapering_frac), dim="time")
    if lowcut is not None and highcut is not None:
        da = xp.filter(
            da, [lowcut, highcut], btype="bandpass", zerophase=False, dim="time"
        )
    elif lowcut is not None and highcut is None:
        da = xp.filter(
            da, lowcut, btype="highpass", zerophase=False, dim="time"
        )
    elif lowcut is None and highcut is not None:
        da = xp.filter(
            da, highcut, btype="lowpass", zerophase=False, dim="time"
        )
    elif lowcut is None and highcut is None:
        pass
    return da


def preprocess_strainrate(da, tapering_frac, lowcut, highcut):
    """
    Function to preprocess displacement data.
    Input:
        da: xdas dadatarray with distance and time
        tapering_frac: fraction of tapering window
        lowcut: low cut frequency in Hz
        highcut: high cut frequency in Hz
        avg_len: window length for m moving average in meters
    Output:
        da: preprocessed displacement dataarray
    """
    da = da - da.mean(dim="distance")
    da = xp.taper(da, window=("tukey", tapering_frac), dim="time")
    if lowcut is not None and highcut is not None:
        da = xp.filter(
            da, [lowcut, highcut], btype="bandpass", zerophase=False, dim="time"
        )
    elif lowcut is not None and highcut is None:
        da = xp.filter(
            da, lowcut, btype="highpass", zerophase=False, dim="time"
        )
    elif lowcut is None and highcut is not None:
        da = xp.filter(
            da, highcut, btype="lowpass", zerophase=False, dim="time"
        )
    elif lowcut is None and highcut is None:
        pass
    return da


def prepare_windows(picks_df, da, pre_win, wlen, noise_wlen, smad_thre):
    """
    Processing steps to isolate starttime of noise and seismic phases.  
    Input:
        picks_df (pd.DataFrame): DataFrame containing picks data.
        da (xdas.DataArray): DataArray with distance and time.
        pre_win (float): Pre-window time in seconds.
        wlen (float): Window length in seconds.
        noise_wlen (float): Noise window length in seconds.
        smad_thre (int): Threshold for SMAD filter.
    Returns:
        starttime (xr.Dataset): Dataset with start times of windows for each phase.
        endtime (xr.Dataset): Dataset with end times of windows for each phase.
        toa (xr.Dataset): Dataset with actual picks for each phase.
    """
    # Generate window timing
    starttime, endtime, toa = window_timing(picks_df, da, pre_win=pre_win, wlen=wlen)
    # Correct window length
    starttime, endtime = correct_pwin(starttime, endtime, wlen)
    # Add noise windows
    starttime, endtime = add_noise(starttime, endtime, da, noise_wlen)
    # Avoid overlaps between phases
    for phases in [["Noise", "P"], ["P", "S"]]:
        if (phases[0] in starttime) & (phases[-1] in starttime):
            starttime, endtime = avoid_overlap(starttime, endtime, phases[0], phases[1])
    # Remove null windows
    for phase in starttime:
        starttime, endtime = avoid_nullwin(starttime, endtime, phase)
    return starttime, endtime, toa


def get_percentile(da_win, perc, dim="time"):
    """
    Function to get the nth percentile amlplitude along given dimension of a dataarray.
    Input:
        da_win: xdas da
        perc: percentile level to calculate
        dim: dimension to calculate the percentile on
    Output:
        max: xdas da with the nth percentile value at each of the other dimension   
        mask: mask containing the index along the given dim at which the percentile is found
    """
    da = np.abs(da_win)
    axis = da.get_axis_num(dim)
    max = np.nanpercentile(da.values, perc, axis=axis)
    mask = np.abs(da - max) == np.nanmin(np.abs(da - max), axis=0)
    mask = mask.to_xarray().drop_duplicates(dim="time")
    return max, mask


def get_win(da, starttime, endtime, phase):
    """
    Function to get the window of a given phase.
    Input:
        da: xdas da
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        phase: phase to get the window for
    Output: 
        da_win: xdas da with values for the selected window and nan values elsewhere
    """
    mask = (da["time"].values.reshape(-1, 1) > starttime[phase].values) & (
        da["time"].values.reshape(-1, 1) < endtime[phase].values
    )
    data = np.full((np.shape(da.values)[0], np.shape(da.values)[-1]), np.nan)
    data[mask] = da.values[mask]
    da_win = da.copy(data=data)
    return da_win


def phases_perc(da, starttime, endtime, perc, dim="time"):
    """
    Function to calculate the nth percentile for each phase.
    Input:
        da: xdas da
        starttime: xarray dataset with start times of windows for each phase
        endtime: xarray dataset with end times of windows for each phase
        perc: percentile level to calculate the peak
        dim: dimension to calculate the percentile on
    Output:
        max: xdas da with the nth percentile value at value each of the other dimension (distance)
        max_mask: mask containing the index along the given dim at which the percentile is found
    """
    max = xr.Dataset()
    max_mask = xr.Dataset()
    for phase in starttime:
        da_win = get_win(da, starttime, endtime, phase)
        max_ph, mask_max_ph = get_percentile(da_win, perc, dim=dim)
        max_ph = np.array(max_ph)
        max_ph = xr.DataArray(max_ph, {"distance": da_win["distance"].values})
        max[phase] = max_ph
        mask_max_ph = xr.DataArray(mask_max_ph, {"time": da_win["time"].values, "distance": da_win["distance"].values})
        max_mask[phase] = mask_max_ph
    return max, max_mask


def phases_snr(max, snr_dict, max_mask, phases=["P", "S"]):
    """
    Function to calculate the signal to noise ratio for each phase.
    Input:
        max: xdas da with the max amplitude estimate along the cable
        snr_dict: dictionary with the SNR threshold for each phase
        max_mask: mask containing the index along the given dim at which the peak is found
        phases: list of phases to calculate the SNR for
    Output:     
        snr: xdas da with the SNR for each phase
        max: xdas da with the max amplitude estimate along the cable
        max_mask: mask containing the index along the given dim at which the peak is found
    """
    snr = xr.Dataset()
    for phase in phases:
        if phase not in max:
            continue
        mask = max[phase] / max["Noise"] > snr_dict[phase]
        snr_ph = max[phase] / max["Noise"]
        snr_ph = xr.DataArray(snr_ph, {"distance": max["distance"].values})
        snr[phase] = snr_ph
        max_mask[phase].values[:, ~mask] = False
    return snr, max, max_mask


def get_hypocentral_distance(event, fiber):
    """
    Calculate the hypocentral distance between an event and a fiber optic cable.
    Input:
        event: pandas df with event information (longitude, latitude, depth in km)
        fiber: pandas df with fiber optic cable information (longitude, latitude, altitude in m)
    Output:
        hyp: xarray dataset with hypocentral distance of each position along the cable
    """
    geod = Geod(ellps="WGS84")
    _, _, epi = geod.inv(
        *np.broadcast_arrays(
            event["longitude"],
            event["latitude"],
            fiber["longitude"],
            fiber["latitude"],
        )
    )
    hyp = np.hypot(epi, event["depth"].values * 1000 + fiber["altitude (m)"].values)
    return xr.DataArray(hyp, {"distance": fiber["distance"].values})


def get_hypodist(event, fiber):
    """
    Calculate the hypocentral distance between an event and a fiber optic cable.
    Input:
        event: pandas df with event information (longitude, latitude, depth in km)
        fiber: pandas df with fiber optic cable information (longitude, latitude, altitude in m)
    Output:
        hyp: xarray dataset with hypocentral distance of each position along the cable
    """
    geod = Geod(ellps="WGS84")
    _, _, epi = geod.inv(
        *np.broadcast_arrays(
            event["longitude"],
            event["latitude"],
            fiber["longitude"].values,
            fiber["latitude"].values,
        )
    )
    hyp = np.hypot(epi, event["depth"] * 1000 + fiber["elevation"].values)
    return xr.DataArray(hyp, {"distance": fiber["distance"].values})


def get_epicentral_distance(event, fiber):
    """
    Calculate the epicentral distance between an event and a fiber optic cable.
    Input:
        event: pandas df with event information (longitude, latitude)
        fiber: pandas df with fiber optic cable information (longitude, latitude)
    Output:
        epi: xarray dataset with epicentral distance of each position along the cable
    """
    geod = Geod(ellps="WGS84")
    _, _, epi = geod.inv(
        *np.broadcast_arrays(
            event["longitude"],
            event["latitude"],
            fiber["longitude"],
            fiber["latitude"],
        )
    )
    return xr.DataArray(epi, {"distance": fiber["distance"].values})


def estimate_ml_ISNET(max, hyp):
    """
    Function to estimate local magnitude from disp DAS data using law from Bobbio et al., 2008.
    Input:
        max: xdas da with the max amplitude estimate along the cable
        snr: xdas da with the signal to noise ratio for each phase along the cable
        hyp: xdas da with hypocentral distance of each position along the cable
        snr_threshold: threshold for signal to noise ratio 
    Output: 
        Mag: xdas da with the local magnitude
    """
    peak = max["S"] * 10 ** -6 # Convert to millimeters
    A = peak * 2800.0
    mag = np.log10(A.values) + 1.79 * np.log10(hyp.values / 1000.0) - 0.58
    return  xr.DataArray(mag, {"distance": max["distance"].values})