#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os,re, sys, math
import obspy
from obspy import read, read_inventory, Trace
import matplotlib.pyplot as plt 
from glob import glob
from scipy import stats
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from math import log10, floor, ceil
import numpy as np
import pandas as pd
from tabulate import tabulate
#for interpolation
#######################################
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Affine2D
#######################
import xlrd
from scipy.signal import detrend
import scipy.signal as signalscp
from scipy.signal import savgol_filter
#################################################
from matplotlib.ticker import FormatStrFormatter, NullFormatter
from scipy.signal import savgol_filter
##############################################
import matplotlib.dates as mdates
#######################################
from matplotlib.colors import Normalize
import matplotlib.cm as cm
############################
import datetime as dtt 
import pandas as pd
#################
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
import matplotlib.dates as mdates
#from pandas.core.common import flatten
#import xarray as xr
#from scipy.signal import convolve2d, find_peaks
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
##################################
from obspy.core import UTCDateTime
from obspy.signal import PPSD
from obspy.imaging.cm import viridis, viridis_white, pqlx
from scipy.signal import welch, csd
#import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from pathlib import Path as Pathdir

AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}
#size for anotation
#xy_coords = (0.71, 0.95)
xy_coords   = (0.85, 0.95)
xy_coords_D = (0.85, 0.8)
#set the fontsize
Figsize= 14
pad    = 20
#define text size
f_size = 14
####for ticks #####
fticks = 14
#############################


def start_endtime(time_serie):
    #minimum time
    tstart     = time_serie.min()
    #tstart     = time_serie[0]
    #Add one day to the minimum time
    tend      = tstart + pd.to_timedelta(1, unit= 'D')
    #tend       = time_serie[-1]
    #Convert to minimum
    tend       = tend.to_numpy()
    return (tstart, tend)
#################################################
#################################################
def read_mseed_only(FILE, xmlFile, comp, Pressure = None):
    #get the stream by reading the mseed file
    st       = read(FILE)
    #get all the channls in the stream
    channels = [tr.stats.channel for tr in st]
    #get the channel you want to extract
    chanel   = channels[0].replace(channels[0][-1], comp)
    #read the inventory
    inv      = read_inventory(xmlFile)
    #prefiltering the signal
    #if(Pressure):
    #    st.remove_response(inventory=inv, pre_filt=pre_filt, output="DEF",water_level=60)
    #else:
    #    #Remove instrument response
    #    st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",water_level=60)
    #sorted the stream
    st.sort(keys=["channel"], reverse=True)
    #print(channels,chanel)
    #extract the trace for the choosing channel
    tr            = st.select(channel=chanel)[0]
    #################################
    return tr







def verical_grid():
    #number of vertical lines
    #num_vertical = 3
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)
    
    return vertical_locations









def FFT_SCP(data, sampling_rate=1.0):
    #Convert seismogram to periods
    freqs, times, Sxx = signalscp.spectrogram(data, fs=sampling_rate, detrend='constant', scaling='density', mode='psd')
    periods = 1.0 / freqs
    #power spectral
    power   = np.mean(Sxx, axis=1)
    return (periods, power)






def TrimTrace(tr, trim_t1, trim_t2):
        #grab the signal and the noise
        #trace_trim    = tr.slice(starttime=start_window[0], endtime=end_window[1]).data
        trace_trim    = tr.slice(starttime=trim_t1, endtime=trim_t2)
        #Calculate RMS for signal and noise
        return trace_trim

def TraceComp(tr):
    if(tr.stats.component=="1"):
        cmpnt = "Z"
    elif(tr.stats.component=="2"):
        cmpnt = "HH2"
    elif(tr.stats.component=="3"):
        cmpnt = "HH3"
    else:
        cmpnt = tr.stats.component
    return cmpnt





## ---------------------------------------------------------
## 1. Load seismic (vertical) and pressure traces
#def remove_coherence(trace_z, trace_p):
#    z = trace_z.data.astype(float)
#    p = trace_p.data.astype(float)
#
#    # Ensure same length
#    n = min(len(z), len(p))
#    z = z[:n]
#    p = p[:n]
#    fs = trace_z.stats.sampling_rate
#    # ---------------------------------------------------------
#    # 2. Compute spectra and transfer function H(f)
#    # ---------------------------------------------------------
#    #nperseg = 4096
#    nperseg = 1024
#    ################
#    noverlap= 2048
#    epsilon = 1e-12
#    coherence_threshold = 0.5
#    #Auto-spectrum of pressure
#    f, Spp = welch(p, fs=fs, nperseg=nperseg)
#    #f, Spp = welch(p, fs=fs, nperseg=nperseg, noverlap=noverlap)
#    #Cross-spectrum between vertical and pressure
#    _, Szp = csd(z, p, fs=fs, nperseg=nperseg)
#    #_, Szp = csd(z, p, fs=fs, nperseg=nperseg, noverlap=noverlap)
#    #Transfer function H(f)
#    H = Szp / Spp
#
##    # Spectra
##    f, Ppp = welch(p, fs, nperseg=nperseg, noverlap=noverlap)
##    _, Pzz = welch(z, fs, nperseg=nperseg, noverlap=noverlap)
##    _, Pzp = csd(z, p, fs, nperseg=nperseg, noverlap=noverlap)
##    H = Pzp / (Ppp + epsilon)
##    # ---------------------------------------------------------
##    # Coherence
##    coh = np.abs(Pzp)**2 / (Ppp * Pzz)
##    # Only trust frequencies with strong coupling
##    H[coh < coherence_threshold] = 0
#    # 3. Apply correction in frequency domain
#    # ---------------------------------------------------------
#    # FFT of signals
#    Zf = np.fft.rfft(z)
#    Pf = np.fft.rfft(p)
#
#    # Interpolate H to FFT frequency grid
#    freq_fft = np.fft.rfftfreq(n, d=1/fs)
#    H_interp = np.interp(freq_fft, f, H.real) + 1j * np.interp(freq_fft, f, H.imag)
#    
#    # Subtract coherent pressure contribution
#    Z_clean_f = Zf - H_interp * Pf
#    
#    #Inverse FFT to time domain
#    z_clean = np.fft.irfft(Z_clean_f, n=n)
#    #---------------------------------------------------------
#    #4. trasform the z_clean into trace
#    #---------------------------------------------------------
#    #station = trace_z.stats.station  
#    #channel  trace_z.stats.channel  
#    #create a new trace
#    #tr = Trace(data = z_clean) # Add metadata 
#    #tr.stats.sampling_rate = fs #Hz 
#    # Create new trace with same metadata 
#    tr       = Trace(data = z_clean)
#    tr.stats = trace_z.stats.copy()
#    ##set the station and the channel
#    #tr.stats.station = station
#    #tr.stats.channel = channel
#    ##get the start time form the initial trace
#    #start = trace_z.stats.starttime 
#    #year  = start.year 
#    #month = start.month 
#    #day   = start.day
#
#    #tr.stats.starttime = UTCDateTime(year, month, day, 0, 0, 0)
#    ##t = np.arange(n) / fs
#    #return(t, z_clean)
#    return(tr)



# 1. Load seismic (vertical) and pressure traces
def remove_coherence(trace_z, trace_p):
    z = trace_z.data.astype(float)
    p = trace_p.data.astype(float)

    # Ensure same length
    n = min(len(z), len(p))
    z = z[:n]
    p = p[:n]
    fs = trace_z.stats.sampling_rate
    # ---------------------------------------------------------
    # 2. Compute spectra and transfer function H(f)
    # ---------------------------------------------------------
    #nperseg = 4096
    nperseg = 1024
    ################
    noverlap= 2048
    epsilon = 1e-12
    coherence_threshold = 0.5
    #Auto-spectrum of pressure
    f, Spp = welch(p, fs=fs, nperseg=nperseg)
    #Cross-spectrum between vertical and pressure
    _, Szp = csd(z, p, fs=fs, nperseg=nperseg)
    #Transfer function H(f)
    H = Szp / Spp

    # 3. Apply correction in frequency domain
    # ---------------------------------------------------------
    # FFT of signals
    Zf = np.fft.rfft(z)
    Pf = np.fft.rfft(p)

    # Interpolate H to FFT frequency grid
    freq_fft = np.fft.rfftfreq(n, d=1/fs)
    H_interp = np.interp(freq_fft, f, H.real) + 1j * np.interp(freq_fft, f, H.imag)
    
    # Subtract coherent pressure contribution
    Z_clean_f = Zf - H_interp * Pf
    
    #Inverse FFT to time domain
    z_clean = np.fft.irfft(Z_clean_f, n=n)
    #---------------------------------------------------------
    #4. trasform the z_clean into trace
    # Create new trace with same metadata 
    tr       = Trace(data = z_clean)
    tr.stats = trace_z.stats.copy()
    return(tr)



#def remove_coherence_sliding_window(trace_z, trace_p, window_len_sec=3600 *2, overlap_frac=0.5):
def remove_coherence_sliding_window(trace_z, trace_p, window_len_sec=3600, overlap_frac=0.5):
    """
    Subtracts coherent pressure noise from the vertical trace using a
    sliding window to account for non-stationary noise.
    """
    z_full = trace_z.data.astype(float)
    p_full = trace_p.data.astype(float)

    # Ensure same length
    n_total = min(len(z_full), len(p_full))
    z_full = z_full[:n_total]
    p_full = p_full[:n_total]

    fs = trace_z.stats.sampling_rate

    # Define window parameters in samples
    win_samples = int(window_len_sec * fs)
    hop_samples = int(win_samples * (1 - overlap_frac))

    # Output array and weight array for overlap-add
    z_clean_full = np.zeros(n_total)
    weights = np.zeros(n_total)

    # Taper to avoid edge effects during overlap-add
    taper = np.hanning(win_samples)

    # Parameters for spectral estimation within the window
    nperseg = min(1024, win_samples // 4)

    # Sliding Window Loop
    for start in range(0, n_total - win_samples + 1, hop_samples):
        end = start + win_samples

        z_win = z_full[start:end]
        p_win = p_full[start:end]

        # 1. Compute local Transfer Function H(f)
        f, Spp = welch(p_win, fs=fs, nperseg=nperseg)
        _, Szp = csd(z_win, p_win, fs=fs, nperseg=nperseg)
        #H = Szp / (Spp + 1e-12) # Added epsilon to prevent div by zero
        H = Szp / Spp  # Added epsilon to prevent div by zero

        # 2. Apply correction in frequency domain for this window
        Zf = np.fft.rfft(z_win)
        Pf = np.fft.rfft(p_win)

        freq_fft = np.fft.rfftfreq(win_samples, d=1/fs)
        H_interp = np.interp(freq_fft, f, H.real) + 1j * np.interp(freq_fft, f, H.imag)

        # Local subtraction
        Z_clean_f = Zf - H_interp * Pf
        z_clean_win = np.fft.irfft(Z_clean_f, n=win_samples)

        # 3. Overlap-Add with tapering
        z_clean_full[start:end] += z_clean_win * taper
        weights[start:end] += taper

    # Handle the end of the trace if it wasn't perfectly divisible by hop_samples
    # Normalize by weights to account for overlapping regions
    idx = weights > 0
    z_clean_full[idx] /= weights[idx]

    # Fill remaining samples (if any) with original data or zeros
    if not np.all(idx):
        z_clean_full[~idx] = z_full[~idx]

    # Create new trace
    tr_clean = Trace(data=z_clean_full)
    tr_clean.stats = trace_z.stats.copy()

    return tr_clean




# 1. Setup paths
# Using glob to grab all mseed files in the directory structure
# You MUST have a stationXML or inventory file for PPSD to work

def calculate_obspy_ppsd(files, inventory_path):
    inv = read_inventory(inventory_path)
    ppsd = None
    
    print(f"Starting PPSD processing of {len(files)} files...")
    
    for i, file_path in enumerate(sorted(files)):
        try:
            st = read(file_path)
            # PPSD needs the trace, but it also needs the statistics from tr.stats
            # We process one trace at a time (e.g., the HH1 channel)
            #tr = st[0]
            
            # Initialize PPSD object on the first successful file
            if ppsd is None:
                # ppsd_length=3600 is standard (1 hour windows)
                ppsd = PPSD(st[0].stats, metadata=inv, ppsd_length=3600)
            
            # Add the current trace data to the PPSD estimate
            ppsd.add(st)
            
            print(f"[{i+1}/{len(files)}] Added: {os.path.basename(file_path)}")
            
            # Clear memory
            del st
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    return ppsd






#Open the configuration file
with open("confreq_coh.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


################
#Grab the window dictionary
WindowingDict     = Fig_params["WindowingDict"]
#Grab the path of the data
#Seimic mseed file
#read the landbased file
mseedFile_Land    = Fig_params['mseedFile_Land']
#Seimic mseed file
Response_File_OBS = Fig_params['Response_File_OBS']
Response_File_Land= Fig_params['Response_File_Land']
#get #Correpondance between OBS componet and landbased station component
OBS_Landbased_comps = Fig_params['OBS_Landbased_comps']
#plot only the ppsd average?
ppsd_average_only   = Fig_params['ppsd_average_only']
#grad the option for y-axis erasizer
Erased_ylabel       = Fig_params['Erased_ylabel']
#grad the option for x-axis erasizer
Erased_xlabel       = Fig_params['Erased_xlabel']
#Erase title?
Erased_title        = Fig_params['Erased_title']
#show noise model level names (NLNM and New High Noise Model NHNM)
show_nm_names       = Fig_params['show_nm_names']
#check if the user want to show the colorbar
show_colorbar       = Fig_params['show_colorbar']
#grab the period limit to visualize as a tuple
period_lim          = tuple(Fig_params['period_lim']) 
#grab the windows
start_1             = Fig_params['start_1']
tend_1              = Fig_params['tend_1']
##Grabe the station type you want to visualte####
Station_Type        = Fig_params['Station_Type'] 
#Set the number of subfigures to Zero

########################################################################
#Set a component you need
component = "HH1" #Vertical compenent of the seismogram
#DEPTH of the seismic station in meter (m) in in the Lake Lucerne
STATIONS_DEPTH = {"MUS03": 59, "MUS06": 71, "MUS07": 69, "MUS05": 49, "MUS04":41, "MUS08": 62, "MUS09": 59, "MUS10": 73}
#open a dictionary
FILEDICT = dict()
path     = Pathdir("SeismicLucerne") 

for item in path.iterdir(): 
    item_new = "%s/mseed"%(item)
    #print(item_new)
    path_years = Pathdir(item_new) 
    for item_year in path_years.iterdir():
        #print(item_year)
        path_months = Pathdir(item_year)
        for item_month in path_months.iterdir():
            #print(item_month)
            path_days = Pathdir(item_month)
            for item_day in path_days.iterdir():
                #print(item_day)
                for File in Pathdir(item_day).glob("*%s*.mseed"%(component)):
                    basename = os.path.basename(str(File))
                    parts    = basename.split("_")
                    station  = parts[1]
                    if station not in FILEDICT:
                        FILEDICT[station] = [] 
                    #else:
                    FILEDICT[station].append(str(File))


#print(FILEDICT)



# Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"
#######################

#####set the figure file name##################################
#filename = "ppsd_%s_%s-%s.pdf"%(sta, dates[0], dates[-1])
##################################################################################################2.Run the calculation

# 1. Initialize the shared figure BEFORE the loop
#import matplotlib.pyplot as plt

# 1. Initialize the shared figure BEFORE the loop

# 1. Initialize the shared figure
fig, main_ax = plt.subplots(figsize=(10, 7))
all_stations_plotted = []

# --- Add Noise Models Manually (Once) ---
periods_nm, nlnm_values = get_nlnm()
periods_nm, nhnm_values = get_nhnm()
main_ax.plot(periods_nm, nlnm_values, color='dimgray', lw=1.5, ls='--', label="NLNM/NHNM", alpha=0.7)
main_ax.plot(periods_nm, nhnm_values, color='dimgray', lw=1.5, ls='--', alpha=0.7)

for sta in FILEDICT:
    FILES_LIST = FILEDICT[sta]
    ppsd_result = calculate_obspy_ppsd(FILES_LIST, Response_File_OBS)

    if ppsd_result:
        # 2. Extract data directly from the PPSD object
        # Note: ObsPy returns (periods, mean_psd_values)
        periods, mean_psd = ppsd_result.get_mean()

        # 3. Plot on the shared axis
        main_ax.plot(periods, mean_psd, lw=2, label=f"Sta: {sta}")
        all_stations_plotted.append(sta)
    else:
        print(f"No data for station {sta}")

# 4. Finalize Formatting
if all_stations_plotted:
    main_ax.set_xscale('log')
    main_ax.set_xlim(period_lim[0], period_lim[1])
    main_ax.set_ylim(-190, -80)

    main_ax.grid(which='both', linestyle=':', color='gray', alpha=0.5)
    main_ax.set_ylabel(r'PSD (dB re $1\ m^2/s^4/Hz$)', fontsize=12)
    main_ax.set_xlabel('Period (s)', fontsize=12)

    main_ax.legend(loc='upper right', fontsize='small', ncol=2)

    if not Erased_title:
        main_ax.set_title("Multi-Station Average PPSD")

    plt.savefig("ppsd_multi_station_average.pdf", dpi=300, bbox_inches="tight")
    #plt.show()
