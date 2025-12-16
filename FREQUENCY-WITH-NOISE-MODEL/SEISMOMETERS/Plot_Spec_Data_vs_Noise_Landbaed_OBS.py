#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os,re, sys, math
import obspy
from obspy import read, read_inventory
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










#Open the configuration file
with open("confreq_noise.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


################
#Grab the window dictionary
WindowingDict     = Fig_params["WindowingDict"]
#Grab the path of the data
#Seimic mseed file
mseedFiles        =   Fig_params['mseedFiles']
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

#fig_space = Fig_params['fig_space']
##Grab the figure size
#fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

########################################################################
#DICTLABEL = {"HH2":None, "HH3":None, "HH1": "Acceleration PSD (m²/s⁴/Hz)", "HDH": "Pressure PSD (Pa²/Hz)", 
#             "Z": None, "N": None, "E": None}
#exit()
for ix,  mseedFile in zip(range(len(mseedFiles)), sorted(mseedFiles)):
    #get the basename of the data file
    basename = os.path.basename(mseedFile)
    #Extract date from the basename using regex
    match    = re.search(r'\d{4}-\d{2}-\d{2}', basename)
    if(match):
        DAYDATE = match.group()
    else:
        print("No date for the file: %s"%(mseedFile))
        exit()
    ###################################################
    cmp_                        = basename.split("-")[0]
    comp                        = basename.split("-")[0].split("_")[-1]
    #get the seismogram deconvolution
    #get the OBS-component
    comp_Land                   = OBS_Landbased_comps[comp]
    if(Station_Type["OBS"]):
        #get the OBS data
        responsexml     = Response_File_OBS
        component       = comp[-1] 
        mseed_file_name = mseedFile
        #start_file_name = "OBS"
    elif(Station_Type["Landbased"]):
        #get the OBS data
        responsexml     = Response_File_Land
        component       = comp_Land 
        mseed_file_name = mseedFile_Land
    ##########################
    inv                         = read_inventory(responsexml)
    stream                      = read(mseed_file_name)
    trace                      = stream.select(component=component)
    #get the date from the trace
    start_date = trace[0].stats.starttime.date
    #transform the start date to a string 
    start_date = str(start_date)
    #get the netwok and station from the trace
    cmpnt = TraceComp(trace[0])
    netwk_stn = "%s.%s.%s"%(trace[0].stats.network, trace[0].stats.station, cmpnt)
    #set start file name
    start_file_name = netwk_stn
    
    ##################################################
    #ppsd                        = PPSD(trace.stats, metadata=inv)
    #########Check windows #########
    if(WindowingDict["Windowing"] == True):
        if(WindowingDict["nwindows"]==1):
            #get the time windows
            start_w1 = UTCDateTime("%sT%s"%(start_date, start_1))
            tend_w1  = UTCDateTime("%sT%s"%(start_date, tend_1))
            #get the trim trace in the window
            tr_w1    = TrimTrace(trace, start_w1, tend_w1)
            ppsd   = PPSD(stats=tr_w1[0].stats, metadata=inv)
            #Add the trace into the ppsd
            ppsd.add(tr_w1)
            ####set the figure file name##################################
            filename = "Fig-psd-%s-%s__%s--%s.pdf"%(start_file_name, start_date, start_1.replace(":","-"), tend_1.replace(":","-"))
        else:
            print("The number of windows should be set to 1")
            exit()
    else:
          ppsd       = PPSD(stats=trace[0].stats, metadata=inv)
          #Add the trace into the ppsd
          ppsd.add(trace)
          ####set the figure file name##################################
          filename = "ppsd-plot-%s-%s.pdf"%(start_file_name, start_date)
######################################
    #fig = ppsd.plot(show=False)
    #periods           = ppsd.period_bin_centers       #Shape: (n_period_bins,)
    #print(periods[:5])
    #exit()
    #set the period limit
    fig = ppsd.plot(show=False, show_noise_models=True, show_coverage=False, show_histogram=True,  show_mode=False,
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
    #main_ax.set_ylabel('PSD (' + r'$\mathrm{10\log_{10}(m^2/s^{4}/Hz}$))', fontsize = Figsize)
    main_ax.set_ylabel('PSD (' + r'$\mathrm{dB \ re \ 1 m^2/s^{4}/Hz}$)', fontsize = Figsize)
    main_ax.set_xlabel('Period (s)', fontsize = Figsize)

    #Customize y-axis tick label font size
    main_ax.tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
    #Label the Y-axis properly
    #main_ax.set_yscale('log')
    main_ax.set_xscale('log')

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
    #get the mean of the periods and mean of the psd and plot them
    periods, mean_psd = ppsd.get_mean()
    #Assuming you already have a PPSD object
    #median_psd = ppsd.get_percentile(percentile=50)
    #main_ax.plot(periods, mean_psd, c='b', lw=10, alpha=0.4)
    main_ax.plot(periods, mean_psd, c='k', lw=1.5, linestyle='-', alpha=0.6)
    main_ax.plot(periods, mean_psd, c='b', lw=1.5, linestyle ='--', alpha=0.8)
    #set the legend on the main_ax
    #legend = main_ax.legend(loc="upper right", fontsize=14, title=netwk_stn, frameon=False,facecolor="blue")
    #Annotate the figure
    main_ax.annotate(
        netwk_stn,                        #Text for the annotation
        xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',         # Use axes fraction for positioning
        fontsize= f_size,                 # Font size for the annotation
        color   = 'k'  ,                  #text color
        alpha   = 0.9 )                   #alpha for the text
        #bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    #set ylimit 
    main_ax.set_ylim(-190, -80)
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
    #Now save the whole figure
    if(WindowingDict["Windowing"] == True):
        fig.axes[0].set_title("Time window: %s—%s (hh:mm) on %s"%(start_1[:-3],tend_1[:-3], start_date),
                          fontsize=9.8, loc="center",fontweight="bold",  color="k")
                          #fontsize=8, loc="center",fontweight="bold",  color="navy")
    else:
        #fig.axes[0].set_title("Date: %s (YYYY:MM:DD)"%(start_date),
        fig.axes[0].set_title("%s"%(start_date),
                          fontsize=14, loc="center",  color="k")
                          #fontsize=14, loc="center",fontweight="bold",  color="k")
    #Erased title from the figure?
    if(Erased_title):
        fig.axes[0].set_title("")
    #Remove background color
    #Remove the colorbar manually
    for ax in fig.axes:
        if hasattr(ax, 'collections') and len(ax.collections) > 0:
            # This is likely the histogram axis; keep it
            continue
        else:
            #Remove any extra axes (like colorbar)
            fig.delaxes(ax)
        #############################
    #############################################
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    #fig.savefig(filename, dpi=300, bbox_inches="tight")
