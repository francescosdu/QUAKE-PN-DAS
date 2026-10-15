import matplotlib.cm as cmx
import matplotlib.colors as colors
#from matplotlib.cm import get_cmap
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from obspy import read, UTCDateTime
import os
from obspy.imaging.spectrogram import spectrogram
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.dates as mdates
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from gamma.utils import association
import pandas as pd
from pyproj import CRS, Transformer
import pandas as pd
import xdas

plt.rcParams['axes.spines.bottom'] = True
plt.rcParams['axes.spines.left'] = True
plt.rcParams['axes.spines.right'] = True
plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.grid.axis'] = 'both'
plt.rcParams['grid.color'] = 'darkgrey'
plt.rcParams['grid.linewidth'] = 0.5
plt.rcParams['axes.axisbelow'] = False
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.ymargin'] = 0.1
# -- Ticks and tick labels --
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['xtick.color'] = 'k'
plt.rcParams['ytick.color'] = 'k'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
# -- Fonts --
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Verdana']
plt.rcParams['text.color'] = 'k'
plt.rcParams['axes.labelcolor'] = 'k'
plt.rcParams['axes.labelsize'] = 12
# -- Saving Options --
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.transparent'] = False
# -- Plot Styles --
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["lines.markersize"] = 3
plt.rcParams["lines.markeredgewidth"] = 1
plt.rcParams["figure.dpi"] = 300

def plot_stations(stream, zoom_left = 0, zoom_right = 120):
    stream.detrend('constant')

    fig, ax = plt.subplots(
        3,
        1,
        figsize=(15, 7),
        sharex=True,
        gridspec_kw={"hspace": 0.05, "height_ratios": [1, 1, 1]},
    )

    for i in range(0,3):
        ax[i].plot(
            stream[i].times(),
            stream[i].data / np.amax(stream[i].data),
            "k",
            label=stream[i].stats.channel,
        )
        ax[i].legend(loc=2)

    ax[1].set_ylabel("Normalised Amplitude")
    ax[2].set_xlabel("Time [s]")
    ax[2].set_xlim(zoom_left, zoom_right)
    plt.show()
    return

def plot_pn_stations(stream, preds, zoom_left = 0, zoom_right = 100):
    color_dict = {"P": "C0", "S": "C1", "Detection": "C2"}
    stream.detrend('constant')

    fig, ax = plt.subplots(
        4,
        1,
        figsize=(15, 7),
        sharex=True,
        gridspec_kw={"hspace": 0.05, "height_ratios": [1, 1, 1, 1]},
    )

    offset = preds[0].stats.starttime - stream[0].stats.starttime
    for pred_trace in preds:
        model, pred_class = pred_trace.stats.channel.split("_")
        if pred_class == "N":
            # Skip noise traces
            continue
        c = color_dict[pred_class]
        ax[3].plot(
            offset + pred_trace.times(), pred_trace.data, label=pred_class, c=c
        )
    ax[3].set_ylabel(model)
    ax[3].legend(loc=2)
    ax[3].set_ylim(-0.01, 1.1)

    for i in range(0,3):
        ax[i].plot(
            stream[i].times(),
            stream[i].data / np.amax(stream[i].data),
            "k",
            label=stream[i].stats.channel,
        )
        ax[i].legend(loc=2)

    ax[1].set_ylabel("Normalised Amplitude")
    ax[3].set_xlabel("Time [s]")
    ax[3].set_xlim(zoom_left, zoom_right)
    plt.show()
    return

def plot_DAS(xr,sat = 10, time_deci=1,space_deci=1,save_name = None):
    mm = 1/25.4
     
    vlim = np.nanmax(np.abs(xr.values))/sat
    print(vlim)
    fig, axes = plt.subplots(
            figsize=(220*mm, 150*mm),
            nrows=1,
            ncols=1,
            sharex=True,
            dpi = 400,
        )

    ax = axes
    xr[::time_deci,::space_deci].plot(
        ax=axes,
        add_colorbar=True,
        add_labels=False,
        norm=mcolors.SymLogNorm(vlim*10**2, vmin=-vlim, vmax=vlim),
        cmap="seismic",
        robust = True,
        cbar_kwargs=dict(
            format=mticker.ScalarFormatter(useMathText=True),
            pad=0,
            shrink=1,
            extend="neither",
            ticks=mticker.MaxNLocator(4),
            label=r"strain rate [$n\varepsilon$/s]")
        )
    ax.set_title('DAS record section', fontsize=20)
    ax.tick_params(bottom=True, left=True, labelbottom=True, labelleft=True)
    ax.set_xlabel("Distance from interrogator [m]")
    ax.set_ylabel("Time [s]")
    ax.set_ylim(xr["time"].values.min(), xr["time"].values.max())
    ax.invert_yaxis()

    ax.xaxis.label.set_size(15)
    ax.yaxis.label.set_size(15)
    ax.title
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    if save_name:
        plt.savefig(save_name +'.png')
    return

def plot_picks(xr,picks,time_deci=1,space_deci=1,sat = 50,size = 0.3, save_name = False, tmin = None, tmax = None):
    mm = 1/25.4
     
    vlim = np.nanmax(xr.values)/sat
    fig, axes = plt.subplots(
            figsize=(220*mm, 150*mm),
            nrows=1,
            ncols=1,
            sharex=True,
            dpi = 400,
        )

    ax = axes
    xr[::time_deci,::space_deci].plot(
        ax=axes,
        add_colorbar=True,
        add_labels=False,
        norm=mcolors.SymLogNorm(vlim*10**2, vmin=-vlim, vmax=vlim),
        cmap="seismic",
        cbar_kwargs=dict(
            format=mticker.ScalarFormatter(useMathText=True),
            pad=0,
            shrink=1,
            extend="neither",
            ticks=mticker.MaxNLocator(4),
            label=r"strain rate [$n\varepsilon$/s]")
        )
    #ax.hlines(np.datetime64(min(picks["time"])),20000,120000)
    ax.set_title(f"DAS record section", fontsize=20)
    ax.tick_params(bottom=True, left=True, labelbottom=True, labelleft=True)
    ax.set_xlabel("Distance from interrogator [m]")
    ax.set_ylabel("Time [s]")
    if tmin:
        ax.set_ylim(pd.to_datetime(tmin), pd.to_datetime(tmax))
    else:
        ax.set_ylim(xr["time"].values.min(), xr["time"].values.max())
    ax.set_xlim(xr["distance"].values.min(), xr["distance"].values.max())
    ax.invert_yaxis()

    ax.xaxis.label.set_size(15)
    ax.yaxis.label.set_size(15)
    ax.title
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    #plot picks
    #picks.insert(5,"index_plot",picks["channel_index"].values * np.mean(np.diff(xr.distance)))
    picks[picks["PHASE"] == "p"].plot(
        ax=axes, x="OFFSET", y="TIME", kind="scatter", color="darkorange", s=size,
        xlabel="Distance from interrogator [m]", ylabel="Time [s]")
    picks[picks["PHASE"] == "s"].plot(
        ax=axes, x="OFFSET", y="TIME", kind="scatter", color="seagreen", s=size, 
        xlabel="Distance from interrogator [m]", ylabel="Time [s]")
    if save_name:
        plt.savefig(save_name +'.png')


def associate(da, phasenet_picks, model_string, gamma_picks, prob_thresh = 0.0, plot = False):
    
    os.makedirs(gamma_picks, exist_ok=True)    
    station_df = pd.read_csv('utils/tito_ch2_geometry_gl_40.csv')
    #correct station name
    station_df["id"] = np.arange(len(station_df)).astype(str)
        

    picks = pd.read_csv(f"{phasenet_picks}/event_{model_string}.txt")
    picks = picks[picks["PROB"] > prob_thresh].reset_index()

    # Projections
    wgs84 = CRS.from_epsg(4326)
    local_crs = CRS.from_epsg(32633)  #Italy 
    transformer = Transformer.from_crs(wgs84, local_crs)

    #
    station_df["x(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[0] / 1e3, axis=1
    )
    station_df["y(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[1] / 1e3, axis=1
    )
    station_df["z(km)"] = -station_df["elevation_m"] / 1e3

    # Gamma
    config = {}
    config["dims"] = ['x(km)', 'y(km)', 'z(km)']
    config["use_dbscan"] = False
    config["use_amplitude"] = False
    config["x(km)"] = (200, 600)  
    config["y(km)"] = (3000, 7000) 
    config["z(km)"] = (-2, 200)
    config["vel"] = {"p": 6.5, "s": 6.5 / 1.75}  # We assume rather high velocities as we expect deeper events

    config["method"] = "BGMM"
    if config["method"] == "BGMM":
        config["oversample_factor"] = 4
    if config["method"] == "GMM":
        config["oversample_factor"] = 1

    config["covariance_prior"] = (100, 100)

    # DBSCAN
    config["bfgs_bounds"] = (
        (config["x(km)"][0] - 1, config["x(km)"][1] + 1),  # x
        (config["y(km)"][0] - 1, config["y(km)"][1] + 1),  # y
        (0, config["z(km)"][1] + 1),  # x
        (None, None),  # t
    )
    config["dbscan_eps"] = 10
    config["dbscan_min_samples"] = 10

        # Filtering
    config["min_picks_per_eq"] = 100
    #config["min_p_picks_per_eq"] = 30
    #config["min_s_picks_per_eq"] = 30
    config["max_sigma11"] = 2.0
    config["max_sigma22"] = 1.0
    config["max_sigma12"] = 1.0


    pick_df = []
    for k in range(0,picks.shape[0]):
        #find correct id in location file
        chann_true = np.where(np.isclose(station_df["distance"],
                                            picks["OFFSET"][k],atol = 2) == True)[0][0]
        pick_df.append({
            "id": str(chann_true),
            "id_old": str(picks["CHANNEL"][k]),
            "timestamp": pd.to_datetime(picks['TIME'][k]),
            "prob": picks['PROB'][k],
            "type": picks['PHASE'][k].lower(),
            "offset": picks['OFFSET'][k]
        })
    pick_df = pd.DataFrame(pick_df)
    #["time"] = pd.to_datetime(df["time"])
    
    #print(pick_df)    
    #print(min(station_df["x(km)"]),max(station_df["x(km)"]))
    #print(min(station_df["y(km)"]),max(station_df["y(km)"]))

    catalogs, assignments = association(pick_df, station_df, config, method=config["method"])

    catalog = pd.DataFrame(catalogs)
    assignments = pd.DataFrame(assignments, columns=["pick_idx", "event_idx", "prob_gamma"])
    #print(catalog)
    if not len(catalog) == 0:
        event_id = catalog["event_index"][np.argmax(catalog["num_picks"])]
        assigments_short = assignments[assignments["event_idx"] == event_id].reset_index()
        pick_ev = pick_df.iloc[assigments_short["pick_idx"]]
        pick_ev.columns = ["ID","ID_OLD","TIME","PROB","PHASE","OFFSET"]
        pick_ev.to_csv(f"{gamma_picks}/event_{model_string}.txt",index = False)
        print(f"Associated events with {catalog['num_p_picks'][np.argmax(catalog['num_picks'])]} P and {catalog['num_s_picks'][np.argmax(catalog['num_picks'])]} S picks")
        if plot:
            plot_picks(da, pick_ev, time_deci = 2, space_deci = 1)
    return 

def associate_2CH(da, phasenet_picks, channel, model_string, gamma_picks, prob_thresh = 0.0, plot = False):
    if not os.path.exists(gamma_picks):
        os.mkdir(gamma_picks)
    
    if channel == 'CH1':
    	station_df = pd.read_csv('tito_ch1_geometry_gl_40.csv')
    else:
        station_df = pd.read_csv('tito_ch2_geometry_gl_40.csv')
    #correct station name
    station_df["id"] = np.arange(len(station_df)).astype(str)
        

    picks = pd.read_csv(f"{phasenet_picks}/event_{model_string}_{channel}.txt")
    picks = picks[picks["PROB"] > prob_thresh].reset_index()

    # Projections
    wgs84 = CRS.from_epsg(4326)
    local_crs = CRS.from_epsg(32633)  #Italy 
    transformer = Transformer.from_crs(wgs84, local_crs)

    #
    station_df["x(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[0] / 1e3, axis=1
    )
    station_df["y(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[1] / 1e3, axis=1
    )
    station_df["z(km)"] = -station_df["elevation_m"] / 1e3

    # Gamma
    config = {}
    config["dims"] = ['x(km)', 'y(km)', 'z(km)']
    config["use_dbscan"] = False
    config["use_amplitude"] = False
    config["x(km)"] = (200, 600)  
    config["y(km)"] = (3000, 7000) 
    config["z(km)"] = (-2, 200)
    config["vel"] = {"p": 6.5, "s": 6.5 / 1.75}  # We assume rather high velocities as we expect deeper events

    config["method"] = "BGMM"
    if config["method"] == "BGMM":
        config["oversample_factor"] = 4
    if config["method"] == "GMM":
        config["oversample_factor"] = 1

    config["covariance_prior"] = (100, 100)

    # DBSCAN
    config["bfgs_bounds"] = (
        (config["x(km)"][0] - 1, config["x(km)"][1] + 1),  # x
        (config["y(km)"][0] - 1, config["y(km)"][1] + 1),  # y
        (0, config["z(km)"][1] + 1),  # x
        (None, None),  # t
    )
    config["dbscan_eps"] = 10
    config["dbscan_min_samples"] = 10

        # Filtering
    config["min_picks_per_eq"] = 100
    #config["min_p_picks_per_eq"] = 30
    #config["min_s_picks_per_eq"] = 30
    config["max_sigma11"] = 2.0
    config["max_sigma22"] = 1.0
    config["max_sigma12"] = 1.0


    pick_df = []
    for k in range(0,picks.shape[0]):
        #find correct id in location file
        chann_true = np.where(np.isclose(station_df["distance"],
                                            picks["OFFSET"][k],atol = 2) == True)[0][0]
        pick_df.append({
            "id": str(chann_true),
            "id_old": str(picks["CHANNEL"][k]),
            "timestamp": pd.to_datetime(picks['TIME'][k]),
            "prob": picks['PROB'][k],
            "type": picks['PHASE'][k].lower(),
            "offset": picks['OFFSET'][k]
        })
    pick_df = pd.DataFrame(pick_df)
    #["time"] = pd.to_datetime(df["time"])
    
    #print(pick_df)    
    #print(min(station_df["x(km)"]),max(station_df["x(km)"]))
    #print(min(station_df["y(km)"]),max(station_df["y(km)"]))

    catalogs, assignments = association(pick_df, station_df, config, method=config["method"])

    catalog = pd.DataFrame(catalogs)
    assignments = pd.DataFrame(assignments, columns=["pick_idx", "event_idx", "prob_gamma"])
    #print(catalog)
    if not len(catalog) == 0:
        event_id = catalog["event_index"][np.argmax(catalog["num_picks"])]
        assigments_short = assignments[assignments["event_idx"] == event_id].reset_index()
        pick_ev = pick_df.iloc[assigments_short["pick_idx"]]
        pick_ev.columns = ["ID","ID_OLD","TIME","PROB","PHASE","OFFSET"]
        pick_ev.to_csv(f"{gamma_picks}/event_{model_string}_{channel}.txt",index = False)
        print(f"Associated events with {catalog['num_p_picks'][np.argmax(catalog['num_picks'])]} P and {catalog['num_s_picks'][np.argmax(catalog['num_picks'])]} S picks")
        if plot:
            plot_picks(da, pick_ev, time_deci = 2, space_deci = 1)
    return 


def plot_compare_auto_manual(manual, auto, traces_DAS, idx_to_plot = 'random', t_left = None, t_right = None):
    #Extract channel to plot
    if idx_to_plot == 'random':
        idx_to_plot = np.random.choice(np.unique(manual["ID"]))

    manual_channel = manual[manual["ID"] == idx_to_plot]
    auto_channel = auto[auto["ID"] == idx_to_plot]

    phases = ['p', 's']
    colors_auto = ['b', 'r']
    colors_man = ['c', 'y']

    fig, ax = plt.subplots(1,1, figsize = (10,4))
    ax.plot(traces_DAS[idx_to_plot - 1].times(), 
            traces_DAS[idx_to_plot - 1].data,
            color = 'k')
    for phase, color in zip(phases,colors_auto):
        if phase in auto_channel["PHASE"].values:
            pick_marker = UTCDateTime(auto_channel["TIME"][auto_channel["PHASE"] == phase].values[0].split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
            ax.vlines(pick_marker, 
                    min(traces_DAS[idx_to_plot - 1].data), 
                    max(traces_DAS[idx_to_plot - 1].data),
                    color = color,
                    linestyles = '--',
                    label = f"{phase.upper()} auto")
    for phase, color in zip(phases,colors_man):
        if phase in manual_channel["phase"].str.lower().values:
            pick_marker = UTCDateTime(manual_channel["time"][manual_channel["phase"] == phase.upper()].values[0].split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
            ax.vlines(pick_marker, 
                    min(traces_DAS[idx_to_plot - 1].data), 
                    max(traces_DAS[idx_to_plot - 1].data),
                    color = color,
                    linestyles = '--',
                    label = f"{phase.upper()} manual")
    edge_left = UTCDateTime(min(auto_channel["TIME"]).split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
    if t_left == None:
        ax.set_xlim(edge_left - 5, edge_left + 10)
    else:
        ax.set_xlim(t_left, t_right)
    ax.set_title(f"Distance: {auto_channel['OFFSET'].values[0]:.1f} m")
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r"DAS amplitude [$n\varepsilon$/s]")
    fig.legend()
    return

def plot_compare_auto_manual_deci(manual, auto, traces_DAS, idx_to_plot = 'random', t_left = None, t_right = None):
    #Extract channel to plot
    if idx_to_plot == 'random':
        idx_to_plot = np.random.choice(np.unique(manual["ID"]))

    manual_channel = manual[manual["ID"] == idx_to_plot]
    auto_channel = auto[auto["ID_OLD"] == idx_to_plot]

    phases = ['p', 's']
    colors_auto = ['b', 'r']
    colors_man = ['c', 'y']

    fig, ax = plt.subplots(1,1, figsize = (10,4))
    ax.plot(traces_DAS[idx_to_plot - 1].times(), 
            traces_DAS[idx_to_plot - 1].data,
            color = 'k')
    for phase, color in zip(phases,colors_auto):
        if phase in auto_channel["PHASE"].values:
            pick_marker = UTCDateTime(auto_channel["TIME"][auto_channel["PHASE"] == phase].values[0].split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
            ax.vlines(pick_marker, 
                    min(traces_DAS[idx_to_plot - 1].data), 
                    max(traces_DAS[idx_to_plot - 1].data),
                    color = color,
                    linestyles = '--',
                    label = f"{phase.upper()} auto")
    for phase, color in zip(phases,colors_man):
        if phase in manual_channel["phase"].str.lower().values:
            pick_marker = UTCDateTime(manual_channel["time"][manual_channel["phase"] == phase.upper()].values[0].split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
            ax.vlines(pick_marker, 
                    min(traces_DAS[idx_to_plot - 1].data), 
                    max(traces_DAS[idx_to_plot - 1].data),
                    color = color,
                    linestyles = '--',
                    label = f"{phase.upper()} manual")
    edge_left = UTCDateTime(min(auto_channel["TIME"]).split('+')[0]) - traces_DAS[idx_to_plot - 1].stats.starttime
    if t_left == None:
        ax.set_xlim(edge_left - 5, edge_left + 10)
    else:
        ax.set_xlim(t_left, t_right)
    ax.set_title(f"Distance: {auto_channel['OFFSET'].values[0]:.1f} m")
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r"DAS amplitude [$n\varepsilon$/s]")
    fig.legend()
    return


def plot_phasenet_das(xr, picks, time_deci=1, space_deci=1, sat = 50, size = 0.3):
    mm = 1/25.4
     
    vlim = np.nanmax(xr.values)/sat
    fig, axes = plt.subplots(
            figsize=(220*mm, 150*mm),
            nrows=1,
            ncols=1,
            sharex=True,
            dpi = 400,
        )

    ax = axes
    xr[::time_deci,::space_deci].plot(
        ax=axes,
        add_colorbar=True,
        add_labels=False,
        norm=mcolors.SymLogNorm(vlim*10**2, vmin=-vlim, vmax=vlim),
        cmap="seismic",
        cbar_kwargs=dict(
            format=mticker.ScalarFormatter(useMathText=True),
            pad=0,
            shrink=1,
            extend="neither",
            ticks=mticker.MaxNLocator(4),
            label=r"strain rate [$n\varepsilon$/s]")
        )
    #ax.hlines(np.datetime64(min(picks["time"])),20000,120000)
    ax.set_title('DAS record section', fontsize=20)
    ax.tick_params(bottom=True, left=True, labelbottom=True, labelleft=True)
    ax.set_xlabel("Distance from interrogator [m]")
    ax.set_ylabel("Time [s]")
    ax.set_ylim(xr["time"].values.min(), xr["time"].values.max())
    ax.invert_yaxis()

    ax.xaxis.label.set_size(15)
    ax.yaxis.label.set_size(15)
    ax.title
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    #plot picks
    picks.insert(5,"index_plot",xr["distance"].values[0] + picks["channel_index"].values * xdas.get_sampling_interval(xr, dim = 'distance'))
    picks[picks["phase_type"] == "P"].plot(
        ax=axes, x="index_plot", y="phase_time", kind="scatter", color="darkorange", s=size,
        xlabel="Distance from interrogator [m]", ylabel="Time [s]")
    picks[picks["phase_type"] == "S"].plot(
        ax=axes, x="index_plot", y="phase_time", kind="scatter", color="seagreen", s=size, 
        xlabel="Distance from interrogator [m]", ylabel="Time [s]")

def associate_DAS(da, picks, gamma_picks, prob_thresh = 0.0, plot = False):
    if not os.path.exists(gamma_picks):
        os.mkdir(gamma_picks)
        
    station_df = pd.read_csv('utils/fiber_coords.csv')
    #correct station name
    station_df["id"] = np.arange(len(station_df)).astype(str)

    picks.columns = ["id", "id_pick","timestamp", "prob", "type", "offset"]
    picks["type"] = picks["type"].str.lower()
    picks = picks[picks["prob"] > prob_thresh].reset_index()

    # Projections
    wgs84 = CRS.from_epsg(4326)
    local_crs = CRS.from_epsg(32633)  #Italy 
    transformer = Transformer.from_crs(wgs84, local_crs)

    #
    station_df["x(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[0] / 1e3, axis=1
    )
    station_df["y(km)"] = station_df.apply(
        lambda x: transformer.transform(x["latitude"], x["longitude"])[1] / 1e3, axis=1
    )
    station_df["z(km)"] = -station_df["altitude (m)"] / 1e3

    # Gamma
    config = {}
    config["dims"] = ['x(km)', 'y(km)', 'z(km)']
    config["use_dbscan"] = False
    config["use_amplitude"] = False
    config["x(km)"] = (200, 600)  
    config["y(km)"] = (3000, 7000) 
    config["z(km)"] = (-2, 200)
    config["vel"] = {"p": 6.5, "s": 6.5 / 1.75}  # We assume rather high velocities as we expect deeper events

    config["method"] = "BGMM"
    if config["method"] == "BGMM":
        config["oversample_factor"] = 4
    if config["method"] == "GMM":
        config["oversample_factor"] = 1

    config["covariance_prior"] = (100, 100)

    # DBSCAN
    config["bfgs_bounds"] = (
        (config["x(km)"][0] - 1, config["x(km)"][1] + 1),  # x
        (config["y(km)"][0] - 1, config["y(km)"][1] + 1),  # y
        (0, config["z(km)"][1] + 1),  # x
        (None, None),  # t
    )
    config["dbscan_eps"] = 10
    config["dbscan_min_samples"] = 10

        # Filtering
    config["min_picks_per_eq"] = 100
    #config["min_p_picks_per_eq"] = 30
    #config["min_s_picks_per_eq"] = 30
    config["max_sigma11"] = 2.0
    config["max_sigma22"] = 1.0
    config["max_sigma12"] = 1.0


    pick_df = []
    for k in range(0,picks.shape[0]):
        #find correct id in location file
        chann_true = np.where(np.isclose(station_df["distance"],
                                            picks["offset"][k],atol = 2) == True)[0][0]
        pick_df.append({
            "id": str(chann_true),
            "id_old": str(picks["id"][k]),
            "timestamp": pd.to_datetime(picks['timestamp'][k]),
            "prob": picks['prob'][k],
            "type": picks['type'][k].lower(),
            "offset": picks['offset'][k]
        })
    pick_df = pd.DataFrame(pick_df)
    #["time"] = pd.to_datetime(df["time"])
    
    #print(pick_df)    
    #print(min(station_df["x(km)"]),max(station_df["x(km)"]))
    #print(min(station_df["y(km)"]),max(station_df["y(km)"]))

    catalogs, assignments = association(pick_df, station_df, config, method=config["method"])

    catalog = pd.DataFrame(catalogs)
    assignments = pd.DataFrame(assignments, columns=["pick_idx", "event_idx", "prob_gamma"])
    #print(catalog)
    if not len(catalog) == 0:
        event_id = catalog["event_index"][np.argmax(catalog["num_picks"])]
        assigments_short = assignments[assignments["event_idx"] == event_id].reset_index()
        pick_ev = pick_df.iloc[assigments_short["pick_idx"]]
        pick_ev.columns = ["ID","ID_OLD","TIME","PROB","PHASE","OFFSET"]
        pick_ev.to_csv(f"{gamma_picks}/event.txt",index = False)
        print(f"Associated events with {catalog['num_p_picks'][np.argmax(catalog['num_picks'])]} P and {catalog['num_s_picks'][np.argmax(catalog['num_picks'])]} S picks")
        if plot:
            plot_picks(da, pick_ev, time_deci = 2, space_deci = 1, sat = 200)
    return 
