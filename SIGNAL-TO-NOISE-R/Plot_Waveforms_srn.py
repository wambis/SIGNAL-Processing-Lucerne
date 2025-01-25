#!/usr/bin/env python
#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
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
#from scipy.stats import signaltonoise
import xlrd
import re
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import xarray as xr
import pandas as pd
from mhkit import dolfyn as dlfn
from datetime import datetime
#################
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



#Read the mseed file
def read_mseed(Date_of_Day, FILE, xmlFile, bandpass, Pressure = None):
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
        #Title         = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
        Title         = "%s  %s  %s"%(network, station, component)
        ##Resampling the data
        tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
        #################################
        #sampling_rate = tr.stats.sampling_rate
        data       = tr.data
        #The sampling rate is 1 second
        date_rng   = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
        #change the date_rng to numpy
        time       = date_rng.to_numpy()
        #tr.times() = time
        #print(tr.times()) 
        #print(time)
        #exit()
    return(time, tr, Title)
    #return(time, data, Title)





def Extract_database(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Perform the Extraction for the correspondig Year and drop all NaN
    dYR  = [df.where(df["YR"]==float(ti.split('-')[0][-2:])).dropna(how='all') for ti in Tlist]
    #Perform the Extraction for the correspondig Month and drop all NaN
    dMTs = [dm.where(dm["MO"]==float(ti.split('-')[1])).dropna(how='all') for dm, ti  in zip(dYR,Tlist)]
    #Perform the Extraction for the correspondig Day and drop all NaN
    dDay = [dd.where(dd["DA"]==float(ti.split('-')[2])).dropna(how='all') for dd, ti  in zip(dMTs, Tlist)]
    #Now let's collect the final Data for the corresponding date for a giving paramter
    #Data_ALL = np.concatenate([ds['Discharge'].to_numpy() for ds in dDay])
    Data_ALL = np.concatenate([ds[param].to_numpy() for ds in dDay])
    #Now let's collect the final time by concatenating the times that has been collected
    Time_ALL = np.concatenate([ds['DateTime'] for ds in dDay])
    #Free Momory by deleting AYR, AMT
    del dYR, dMTs
    #Return the  value corresponding to the date of your choice
    return (Time_ALL, Data_ALL)

def Add_Column(ds):
    ds['DateTime'] = pd.to_datetime(ds[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
    #Return the dataframe
    return ds

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




##Get the parameters
fsize = 14
pad   = 10

##############################################
#Define function
def Plot_fig(ax, ix, time, data, Title, waveform= None):
    if(waveform):
        #ylabel = 'ACC (µm/s²)'
        #ylabel = 'Normalized Amp (A/Amax)'
        #ylabel = r"$\frac{A}{A_{\mathrm{max}}}$"
        ylabel = r"$A/A_{\mathrm{max}}$"
        #ylabel = r"$A/A_max$"
    else:
        ylabel = 'Pressure (Pa)'
    color      = "k"
    ylabel     = re.sub('', "", ylabel)
    #ax.plot(time, data, kwargs)
    network, station, component = Title.split("_")
    param      = "%s--%s"%(station, component)
    #get the start and the endtime
    tstart, tend =  start_endtime(time)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0, label = param.title())
    #ax.plot(time, data, lw=0.4, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label =param)
    ax.plot(time, data/max(data), lw=0.3, linestyle ="-", color = color, alpha =0.9, label =param)

    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #get the coordinate of the point, where to write the text
    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    #ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=(0.0089, 0.912),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.87),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
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
    ax.set_ylim(-1.0, 1.0)
#    ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])










def SignalTonoiseFull(data):
    #Compute the signal-to-noise ratio
    data   = data/max(data)
    #snr = signaltonoise(data)
    mean    = np.mean(data)
    std_dev = np.std(data)
    snr     = mean / std_dev
    return abs(snr)


def SignalTonoiseR(tr, signal_window, noise_window):
        #grab the signal and the noise
        signal    = tr.slice(starttime=signal_window[0], endtime=signal_window[1]).data
        noise     = tr.slice(starttime=noise_window[0], endtime=noise_window[1]).data
        #Calculate RMS for signal and noise
        signal_rms= np.sqrt(np.mean(signal**2))
        noise_rms = np.sqrt(np.mean(noise**2))
        #Calculate SNR
        snr = signal_rms / noise_rms
        return snr


def plot_FFT(ax, tr):

   ##Determine the frequencies
    sampling_rate = tr.stats.sampling_rate
    npts          = tr.stats.npts

    #Compute the FFT
    FFT_w = np.fft.fft(tr.data)
    #Spectrum
    spectrum  = ( np.abs(FFT_w) ** 2 * (npts/sampling_rate) )

    #compute the frequencies
    freqs = np.linspace(0, sampling_rate, npts, endpoint= False)



    #Plot the FF transform 
    ax.plot(freqs *1000, spectrum/max(spectrum), c='k', lw=.5)

    #Set Title
    ax.set_title('Fourrier-Transform', fontsize=12, loc ='center', pad =.2)
    #speed up the plot
    plt.style.use('fast')
    #ax.axis('tight')
    ax.set(xlabel='frequency (mHz)', ylabel='PSD (W/Hz)')
    fmax = 5.0
    ax.set_xlim([0, fmax])
    #ax.set_xlim([994, 1000])






#Open the configuration file
with open("confs.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFiles      =   Fig_params['mseedFiles']
STARTDATE       =   Fig_params['STARTDATE']
#Seimic mseed file
Response_File  = Fig_params['Response_File']
#Grab the path of the data
#############################################
#Check the user want to make a plot bigger?
#plot_bigger    = Fig_params['plot_bigger']
# Plot the seismic Waveform 
waveform       = Fig_params['waveform']
pressure_waveform= Fig_params['pressure_waveform']
#Grab the bandpass
bandpass       = Fig_params['bandpass']
#Calculate the signal_to_noise_ratio
signal_to_noise_r = Fig_params['signal_to_noise_r']
#get the window

startime_noise  = Fig_params['startime_noise']
endtime_noise   = Fig_params['endtime_noise']
#set the start and the endtime for signal
startime_signal = Fig_params['startime_signal']
endtime_signal  = Fig_params['endtime_signal']

#Set the number of subfigures to Zero
#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs     = len(mseedFiles)
#Create the figure
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])


#Create the figure
#fig, axs= plt.subplots(nfigs, 1, sharex = False, figsize=(12,10))

fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#Set the window 
# Define signal and noise windows (replace with your own times)
noise_window = (UTCDateTime("%sT%s"%(STARTDATE, startime_noise)), UTCDateTime("%sT%s"%(STARTDATE, endtime_noise)))
###signal_window = (UTCDateTime("2025-01-24T00:00:00"), UTCDateTime("2025-01-24T00:01:00"))
signal_window= (UTCDateTime("%sT%s"%(STARTDATE, startime_signal)), UTCDateTime("%sT%s"%(STARTDATE,endtime_signal)))
#######Plot the discharge ##############
#plot_discharge(ax1,hours_of_day, data_day, starttime, Remove_frame)

#Diction for signal-to-noise ratio
DictSNR = dict()
for ix,  mseedFile in zip(range(len(mseedFiles)), mseedFiles):
    #plot_bigger = False
    basename = os.path.basename(mseedFile)
    cmp_     = basename.split("-")[0]
    #get the seismogram from mseedFile
    time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass, Pressure = False)
    #network, station, component = Title.split()
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], ix, time, tr.data, cmp_, waveform=True)
    #normalized the data
    tr.data = tr.data/max(tr.data)
    snr = SignalTonoiseR(tr, signal_window, noise_window)
    DictSNR[ix] = float("%.2f"%(snr))
    print('%s : %s'%(mseedFile, str(snr)))


#Set the axis-label for the last axis
#format='%y %m %d %H %M %S'
axs[nfigs - 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#axs[nfigs -1].xaxis.set_major_formatter(DateFormatter('%H:%M'))
basename_ = os.path.basename(mseedFile) 
basename = basename_.split(".mseed")[0]

    #Set the font size of yticks
startime_noise  = Fig_params['startime_noise']
endtime_noise   = Fig_params['endtime_noise']
#set the start and the endtime for signal
startime_signal = Fig_params['startime_signal']
endtime_signal  = Fig_params['endtime_signal']

#wind1_t1 = pd.Timestamp('2023-10-23 07:12')
###############Noise window ############################
wind1_t1 = pd.Timestamp('%s %s'%(STARTDATE,startime_noise))
wind1_t2 = pd.Timestamp('%s %s'%(STARTDATE, endtime_noise))
###############Signal window ############################
wind2_t1 = pd.Timestamp('%s %s'%(STARTDATE, startime_signal))
wind2_t2 = pd.Timestamp('%s %s'%(STARTDATE, endtime_signal))
##############Time to write the SNR %%%%%%%%%%%%%%%%%%
time_write_2 = pd.Timestamp('%s 04:48'%(STARTDATE))
for ii in range(nfigs):
    #noise's window
    axs[ii].axvline(wind1_t1, color='k', linestyle='--', alpha=0.4)
    axs[ii].axvline(wind1_t2, color='k', linestyle='--', alpha=0.4)
    #signal's window
    axs[ii].axvline(wind2_t1, color='r', linestyle='--', alpha=0.4)
    axs[ii].axvline(wind2_t2, color='r', linestyle='--', alpha=0.4)
    #make a window between the time
    axs[ii].axvspan(wind1_t1, wind1_t2, color='b', alpha=0.3)
    axs[ii].axvspan(wind2_t1, wind2_t2, color='r', alpha=0.3)
    ####Write on the axis #########
    #axs[ii].text(time_write_2, 0.45, 'SNR=%d'%(DictSNR[ii]), fontsize=15, color='k')
    axs[ii].text(time_write_2, 0.45, 'SNR=%.1f'%(DictSNR[ii]), fontsize=15, color='k')

#Write text on the x-axis at specific_time
#axs[nfigs-1].text(time_write_1, 2.0, 'Wind onset', rotation=0, color='k')
#axs[nfigs-1].text(time_write_2, 2.0, 'Seismic signal onset', rotation=0, color='k')

#### Write on the x-axis
axs[nfigs -1].set_xlabel('Time (hour:minute) on %s'%(STARTDATE), labelpad= pad, fontsize = fsize)
figname   = "Waveforms_%s_of_%s.png"%(basename, STARTDATE)
#figname   = "Waveforms_%s_of_%s.pdf"%(basename, STARTDATE)
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
#fig.savefig(figname)
