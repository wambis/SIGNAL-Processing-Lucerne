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

from pathlib import Path as Pathdir

AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}
#size for anotation
#xy_coords = (0.0081, 0.90)
xy_coords = (0.71, 0.95)
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







def FFT(data, sampling_rate=1):
    """
    Perform FFT on velocity data and plot frequency spectrum.
    Parameters:
        velocity_data (array-like): Velocity data for analysis.
        sampling_rate (float): Sampling rate in Hz.
    Returns:
        freqs (numpy.ndarray): Frequencies of the spectrum.
        power (numpy.ndarray): Power of the frequency spectrum.
    """
    #set the sampling_rate to
    #sampling_rate = 1.0
    #Detrend the data to remove linear trends
    mean = np.mean(data) 
    data = data - mean
    data = detrend(data)
    #Aplly a window function to reduce spactral leakage
    #window = np.hanning(len(data))
    #data   = data * window
    #Perform FFT
    n    = len(data)
    fft_result = np.fft.fft(data)
    freqs = np.fft.fftfreq(n, d=1/sampling_rate)
    power = np.abs(fft_result)**2

    #Keep only positive frequencies
    positive_freqs = freqs[:n // 2]
    positive_power = power[:n // 2]
    ##sorted find the index of the highest peak
    #p_index = np.argmin(positive_power)
    ##Arrange the data such as that the highest peak corresponds to the maximum energy
    #sorted_index =np.argsort(positive_power)
    #freq_sorted  = positive_freqs[sorted_index]
    #power_sorted  = positive_power[sorted_index]

    #positive_freqs = freq_sorted
    #positive_power = power_sorted
    positive_freqs = 1./positive_freqs 
    #return value
    return (positive_freqs, positive_power)
    #return (positive_freqs, positive_power/max(positive_power))
    #return (positive_freqs *1000, positive_power)



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

#def calculate_obspy_ppsd(files, inventory_path):
#    inv = read_inventory(inventory_path)
#    ppsd = None
#    
#    print(f"Starting PPSD processing of {len(files)} files...")
#    
#    for i, file_path in enumerate(sorted(files)):
#        try:
#            st = read(file_path)
#            # PPSD needs the trace, but it also needs the statistics from tr.stats
#            # We process one trace at a time (e.g., the HH1 channel)
#            #tr = st[0]
#            #manually remove the response for hydrophone pressure stream that will be add in the ppsd
#            st.remove_response(inventory=inv, output="DEF", pre_filt=(0.005, 0.006, 45, 50), zero_mean=True, water_level=60)
#            # Initialize PPSD object on the first successful file
#            if ppsd is None:
#                ##################################################
#                #tr.remove_response(inventory=inv, output="DEF", pre_filt=(0.001, 0.002, 45, 50), zero_mean=True, water_level=60)
#                # ppsd_length=3600 is standard (1 hour windows)
#                ppsd = PPSD(stats=st[0].stats, metadata=inv, ppsd_length=3600, special_handling="hydrophone")
#                #ppsd = PPSD(stats=tr.stats, metadata=inv, special_handling="hydrophone")
#            
#            # Add the current trace data to the PPSD estimate
#            ppsd.add(st)
#            
#            print(f"[{i+1}/{len(files)}] Added: {os.path.basename(file_path)}")
#            
#            # Clear memory
#            del st
#            
#        except Exception as e:
#            print(f"Error processing {file_path}: {e}")
#
#    return ppsd



def calculate_obspy_ppsd(files, inventory_path):

    inv         = read_inventory(inventory_path)

    first_trace = read(files[0])[0]

    ppsd        = PPSD(stats=first_trace.stats, metadata=inv, ppsd_length=3600, special_handling="hydrophone")

    for i, file_path in enumerate(files):

        try:
            tr = read(file_path)[0]

            tr.remove_response(
                inventory=inv,
                output="DEF",
                pre_filt=(0.005,0.006,45,50),
                zero_mean=True,
                water_level=60)

            ppsd.add(tr)

            if i % 50 == 0:
                print(f"{i}/{len(files)} processed")

        except Exception as e:
            print(f"Error {file_path}: {e}")

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
#grad the option for y-axis erasizer
Erased_ylabel       = Fig_params['Erased_ylabel']
#grad the option for x-axis erasizer
Erased_xlabel       = Fig_params['Erased_xlabel']
#Erase title?
Erased_title         = Fig_params['Erased_title']
#show noise model level names (NLNM and New High Noise Model NHNM)
show_nm_names       = Fig_params['show_nm_names']
#check if the user want to show the colorbar
show_colorbar       = Fig_params['show_colorbar']
#grab the period limit to visualize as a tuple
period_lim          =tuple(Fig_params['period_lim']) 
#grab the windows
start_1             = Fig_params['start_1']
tend_1              = Fig_params['tend_1']
##Grabe the station type you want to visualte####
Station_Type        = Fig_params['Station_Type'] 
#Set the number of subfigures to Zero

########################################################################
#             "Z": None, "N": None, "E": None}
#exit()
FILES_LIST = []

path = Pathdir("SeismicLucerne") 

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
                for File in Pathdir(item_day).glob("*HDH*.mseed"):
                    #Files = glob(item_day+"/*.mseed")
                    #print(File)
                    ##############
                    FILES_LIST.append(File)




# Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"
#######################
dates = [re.search(date_pattern, str(path)).group() for path in FILES_LIST if path and re.search(date_pattern, str(path))]
######################
dates = sorted(set(dates))

#####set the figure file name##################################
filename = "ppsd-Pressure-%s-%s.pdf"%(dates[0], dates[-1])
#2.Run the calculation
ppsd_result = calculate_obspy_ppsd(FILES_LIST, Response_File_OBS)
#3. Plotting
#ObsPy's PPSD has a built-in plotting method that is much more informative
if ppsd_result:
    ########################## 
    fig = ppsd_result.plot(show=False, show_noise_models=True, show_coverage=False, show_histogram=True,  show_mode=False,
                      show_mean=False, show_earthquakes=None,grid=True, cumulative=False,cumulative_number_of_colors=200,
                      cmap=pqlx, period_lim=period_lim)
    #get all axis
    axes    = fig.axes
    main_ax = axes[0]
    # Customize major gridlines
    main_ax.grid(which='both', linestyle='--', color='gray', linewidth=0.5)
    #set the main_ax limit
    main_ax.set_xlim(period_lim[0],period_lim[1])
    #Coustumize main_ax label
    #main_ax.set_ylabel(r'$\mathrm{dB \ re \ 1 m^2/s^{4}/Hz}$', fontsize = Figsize)
    main_ax.set_ylabel('PSD (' + r'$\mathrm{dB \ re \ 1 Pa^{2}/Hz}$)', fontsize = Figsize)
    main_ax.set_xlabel('Period (s)', fontsize = Figsize)
    
    #Customize y-axis tick label font size
    main_ax.tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
    #Label the Y-axis properly
    #main_ax.set_yscale('log')
    main_ax.set_xscale('log')
    
    
    #Customize y-axis tick label font size
    main_ax.tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
    #Label the Y-axis properly
    #main_ax.set_yscale('log')
    main_ax.set_xscale('log')
    
    #Erased title from the figure?
    if(Erased_title):
        fig.axes[0].set_title("")
    else:
        fig.axes[0].set_title("%s--%s"%(dates[0],dates[-1]),
                          fontsize=14, loc="center",  color="k")
    #erased the ylabel if needed
    if(Erased_ylabel):
        #disable ticks on  y-axis
        #main_ax.yaxis.set_visible(False)
        #main_ax.set_yticks([])        # remove tick marks
        main_ax.set_yticklabels([])   # remove tick labels
        main_ax.set_ylabel("")        # remove the axis label
    if(Erased_xlabel):
        #disable ticks on  x-axis
        #main_ax.set_xticks([])        # remove tick marks
        main_ax.set_xticklabels([])   # remove tick labels
        main_ax.set_xlabel("")        # remove the axis label
    
    ################################

    #Check if the user want the see the colorbar
    if(show_colorbar):
        cbar = fig.axes[-1] 
        cbar.set_ylabel('Probability (%)')
        #ticks = cbar.get_yticks()
        #cbar.set_yticklabels([f"{tick/100:.2f}" for tick in ticks])
        #cbar.set_ylabel("Probability")
        for ax in axes[1: -1]:
            ax.set_visible(False)
            fig.delaxes(ax)
    else:
        #if you're plotting for several day do not delete it
        for ax in axes[1:]:
            ax.set_visible(False)
            fig.delaxes(ax)
    #set a desired color for New Low Noise Model (NLNM) and New High Noise Model (NHNM) 
    color_nm = "dimgray"
    fntsize  = 12
    #The noise models are always the last two lines added, we can then do:
    lines = main_ax.get_lines()
    lines[-2].set_color(color_nm)   # Likely NLNM
    lines[-1].set_color(color_nm)  # Likely NHNM
    #ydata = lines[-2].get_ydata()
    #set the linewidth
    lines[-2].set_linewidth(1.5)         #Your desired linewidth
    lines[-2].set_linestyle("-")         #Solid line
    lines[-2].set_alpha(1.0)             #set the alpha
    #
    lines[-1].set_linewidth(1.5)         #Your desired linewidth
    lines[-1].set_linestyle("-")         #Solid line
    lines[-1].set_alpha(1.0)             #set the alpha
    #Check if the user want to add NLNM and NHNM  at a specific location (x = period, y = PSD)
    if(show_nm_names):
        main_ax.text(0.25, -165, "NLNM", color=color_nm, fontsize=fntsize)
        main_ax.text(0.25, -100, "NHNM", color=color_nm, fontsize=fntsize)
    #set ylimit 
    #main_ax.set_ylim(-160, -60)
    #get the mean of the periods and mean of the psd and plot them
    periods, mean_psd = ppsd_result.get_mean()
    #Assuming you already have a PPSD object
    #median_psd = ppsd_result.get_percentile(percentile=50)
    #main_ax.plot(periods, mean_psd, c='b', lw=10, alpha=0.4)
else:
    print("No data was processed.")


###############################################
fig.savefig(filename, dpi=300, bbox_inches="tight")

#######################################
