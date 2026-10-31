##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
from scipy.signal import welch
from scipy import signal
import os, re, json, yaml
from yaml.loader import SafeLoader
import xlrd, sys, math
import gc
from obspy import read, read_inventory, Stream
import matplotlib.ticker as mticker
from math import log10, floor, ceil
import numpy as np
import pandas as pd
from tabulate import tabulate
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
import matplotlib.dates as mdates
from mhkit.dolfyn.adp import api
from pandas.core.common import flatten
##################################
from obspy.core import UTCDateTime


AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            5: 'f', 6: 'f', 7:'h', 8:'g', 9: 'j',10:'k',
            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}


xy_coords = (0.0065, 0.92) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.95) # Coordinates (relative to axes) for the annotation


#Set parameters
#xy_coords   = (0.85, 0.95)
#xy_coords_D = (0.55, 0.85)
#set the fontsize
Figsize= 14
pad    = 10
#define text size
f_size = 14
##Get the parameters
fsize = 14
####for ticks #####
fticks = 14


def compute_coherence(tr1, tr2, fs):
    # ---------------------------------------------------------
    # 2. Compute PSD using Welch method (Nakamura requires PSD)
    # ---------------------------------------------------------
    #nperseg = 256  # adjust depending on your sampling rate
    #nperseg = 250000  # adjust depending on your sampling rate
    T_max = 1000  #second
    #Compute the nperseg as function sampling rate and the maximum periods
    nperseg = int(fs * T_max)
    
    # ---------------------------------------------------------
    # 3. Compute Horizontal Spectrum (root-mean-square)
    f, Cxy = signal.coherence(tr1, tr2, fs, nperseg=nperseg)
    #f, Cxy = signal.coherence(trZ, trP, fs, nperseg=nperseg, nfft=3600)
    ######################################################
    # Remove f = 0 (avoid division by zero)
    mask = f > 0
    f = f[mask]
    Cxy = Cxy[mask]
    
    # Convert frequency to period
    T = 1.0 / f
    
    # Select period range between 1e-3 and 1000 seconds
    mask_period = (T >= 1e-3) & (T <= 1000)
    T = T[mask_period]
    Cxy = Cxy[mask_period]
    return (T, Cxy)



#Read the mseed file
def read_mseed(Date_of_Day, FILE, xmlFile, bandpass, Pressure=None):
    #get the stream by reading the mseed file
    st      = read(FILE)
    #read the inventory
    inv     = read_inventory(xmlFile)
    #prefiltering the signal
    #Very Good match with the discharge of the day with this event  XJ_MUS01_HH2-2018-12-24.mseed
    dT      = '%s s'%(bandpass)
    #Periods List
    pds     = bandpass.split("-")
    f1      = 1./float(pds[-1]); f2 = 1./float(pds[-2]) ;
    ####################################################
    f3      = 1./float(pds[1]);  f4 = 1./float(pds[0])
    #f1 = 0.001; f2 = 0.00125 ; f3 = 0.002; f4 = 0.0025
    pre_filt = [f1, f2, f3, f4]
    if(Pressure):
        st.remove_response(inventory=inv, pre_filt=pre_filt, output="DEF",water_level=60)
    else:
        #Remove instrument response
        st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",water_level=60)
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",water_level=60)

           #get the trace from the stram
    for tr in st:
        #get the parameters of the trace
        network       = tr.stats.network
        station       = tr.stats.station
        component     = tr.stats.component
        channel       = tr.stats.channel
        stime         = str(tr.stats.starttime)
        endtime       = str(tr.stats.endtime)
        #remove mean and trend on the trace
        tr.detrend(type='demean')
        if(component=="1"):
            component ="Z"
        elif(component=="2"):
            #component ="HH2"
            component ="X"
        elif(component=="3"):
            #component ="HH3"
            component ="Y"
        else:
            component ="HD"
        #Title         = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
        #Title         = "%s  %s  %s"%(network, station, component)
        Title         = "%s.%s.%s"%(network, station, component)
        ##Resampling the data
        tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
        #tr.resample(sampling_rate = 2.0, window = 'hann', no_filter = True)
        #################################
        #The sampling rate is 1 second
        date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= tr.data.size)
        #date_rng  = pd.date_range(start= Date_of_Day, freq='0.5s', periods= tr.data.size)
       #wamba
        #Create datetime array
        ## Convert to pd.date_range (optional, if you want a fixed frequency index)
#        date_rng = pd.date_range(start=Date_of_Day, periods=tr.data.size, freq=f"{tr.stats.delta}S")
        #change the date_rng to numpy
        time      = date_rng.to_numpy()
        #print(tr.data.size, len(time))
        #exit()
    #return(time, data, Title)
    return(time, tr, Title)


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


def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations




#Define function
def Plot_fig(ax, ix, time, data, Title, waveform= None):
    if(waveform):
        ylabel = 'ACC (µm/s²)'
        #ylabel = 'Normalized Amp (A/Amax)'
        #ylabel = r"$\frac{A}{A_{\mathrm{max}}}$"
        #ylabel = r"$A/A_{\mathrm{max}}$"
        #ylabel = r"$A/A_max$"
    else:
        ylabel = 'Pressure (Pa)'
    color      = "k"
    ylabel     = re.sub('', "", ylabel)
    #ax.plot(time, data, kwargs)
    network, station, component = Title.split(".")
    param      = "%s--%s"%(station, component)
    #get the start and the endtime
    tstart, tend =  start_endtime(time)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label =param)
    #ax.plot(time, data/max(data), lw=0.35, linestyle ="-", color = color, alpha =0.9, label =param)
    ##########
    ax.plot(time, data, lw=0.35, linestyle ="-", color = color, alpha =0.9, label =param)
    #ax.plot(time, data, lw=0.5, linestyle ="-", color = color, alpha =0.8, label =param)
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        #get the coordinate of the point, where to write the text
    #Write the text on the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy= xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=16,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #Add the legend
    ax.legend(loc="upper right")
    #ax.legend(loc="upper left", fontsize='large',  bbox_to_anchor=(0, 1), frameon=True, shadow=False)
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = fsize)
    #ax.yticks(fontsize = fsize)
    ax.tick_params(axis='y', labelsize=fsize)
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
#    ax.set_ylim(-1.0, 1.0)
#    ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])





#config_file = "configs.yaml"
#Open the configuration file
with open("confcom.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)

# ---------------------------------------------------------
# 1. Load your three-component seismic data



#Seimic mseed file
mseedFiles          =   Fig_params['mseedFiles']
###########################################
#Seimic mseed file
Response_File  = Fig_params['Response_File']
#get the filter band
bandpass       = Fig_params['bandpass']
# ---------------------------------------------------------
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

##set params values
#xlabel     = 'Period (s)'
##plt.ylabel('Coherence')
#ylabel_sis = 'ACC (µm/s²)'
#ylabel_p   = 'Pressure (Pa)'
####################
nfigs      = len(mseedFiles)
######################
##################Create a figure#######
fig, axs   = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
# Replace with your own file or stream

#loop over the mseedFiles List
for ix,  mseedFile in zip(range(len(mseedFiles)), sorted(mseedFiles)):
    #plot_bigger = False
    basename = os.path.basename(mseedFile)
    station  = basename.split("_")[1]
    #Extract date from the basename using regex
    match    = re.search(r'\d{4}-\d{2}-\d{2}', basename)
    if(match):
        DAYDATE = match.group()
    else:
        print("No date for the file: %s"%(mseedFile))
        exit()
    ###################################################
    cmp_     = basename.split("-")[0]
    #set the parameter
    #get the seismogram from mseedFile
    #time, tr, Title       = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure)
    #get the seismogram from mseedFile
    if("HDH" in cmp_):
        time, tr, Title   = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = True)
        #split
        network, station, component = Title.split(".")
        #Plot the figure
        Plot_fig(axs[ix], ix, time, tr.data, Title)
        #axs[ix].set_ylim(-5, 5) #for 2023-09-22
    else:
        #time, data, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = False)
        time, tr, Title   = read_mseed(DAYDATE, mseedFile, Response_File, bandpass)
        ####################
        network, station, component = Title.split(".")
        data = tr.data * 1e+6
        Plot_fig(axs[ix], ix, time, data, Title, waveform=True)
        ###################
        #if("HH1" in cmp_):
        #axs[ix].set_ylim(-1.0, 1.0) #Good for 2023-08-28

figname = "Waveform_comp_%s.%s.pdf"%(network ,station)

#### Write on the x-axis
axs[nfigs - 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[nfigs -1].set_xlabel('Time (hour:minute) on %s'%(DAYDATE), labelpad= pad, fontsize = fsize)
##########################################################
#Set the font size of yticks
#plt.yticks(fontsize=fsize)
# Set the font size of xticks
plt.xticks(fontsize=fsize)
#Make a rotation of xticks
#plt.xticks(rotation=45)
#Space between the subplots
#plt.subplots_adjust(hspace = 0.08)
plt.subplots_adjust(hspace = fig_space)
#Align ylabel of the figure
fig.align_ylabels()
#Save the figure
#figname = "TEST-PLOTS.png"
fig.savefig(figname, bbox_inches = 'tight')



##set ylimit
##write the text on the figure
##Remove the ticks marks (small vertical lines) on the x-axis
##ax.tick_params(axis="both", length = 0, color= "white", width = 0)
#ax.tick_params(axis='both', labelsize=Figsize)
#ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
#ax.set_xlabel(xlabel, labelpad = pad, fontsize = Figsize)
##save the figure
#
#figname = "Coherence-%s-%s.pdf"%(station, dates[0])
##Save the figure
#plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
