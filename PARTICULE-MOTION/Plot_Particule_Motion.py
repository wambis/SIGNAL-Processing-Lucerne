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
from scipy.interpolate import InterpolatedUnivariateSpline
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


##Get the parameters
fsize = 14
#############
ftext = 18
#fsize = 11
pad   = 10
################################
#xy_coords = (0.0065, 0.85) #Coordinates (relative to axes) for the annotation
#xy_coords = (0.0099, 0.83) #Coordinates (relative to axes) for the annotation
xy_coords = (0.013, 0.83) #Coordinates (relative to axes) for the annotation

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
        data          = tr.data
        #time          = tr.times()
        #The sampling rate is 1 second
        date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
        #change the date_rng to numpy
        time      = date_rng.to_numpy()
        #time, data    = InterSpline(Date_of_Day, tr.times(), tr.data)
        #return(time, tr, Title)
    return(time, data, Title)





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
    network, station, component = Title.split("_")
    param      = "%s--%s"%(station, component)
    #get the start and the endtime
    tstart, tend =  start_endtime(time)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label =param)
    #ax.plot(time, data/max(data), lw=0.35, linestyle ="-", color = color, alpha =0.9, label =param)
    ax.plot(time, data, lw=0.35, linestyle ="-", color = color, alpha =0.9, label =param)

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
        fontsize=16                 # Font size for the annotation
        #bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        #bbox=dict(boxstyle='square', facecolor='white', edgecolor='white', alpha=0.001)  # Box style
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











def Plot_fig2D(ax,fig_object, data2D, **PARAMSDICT):
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    vmind = np.min(data2D, axis = 1) 
    vmaxd = np.max(data2D, axis = 1)
    #Make a 2D plot
    #im = ax.imshow(data2D, cmap = color, vmin = 0. , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
    im = ax.imshow(data2D, cmap = color, vmin = min(vmind) , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
                   interpolation ='spline16')
    #ax.imshow(data2D, cmap = color, origin = 'lower', aspect = 'auto')
    # add the colorbar using the figure's method,
    #telling which mappable we're talking about and
    #which axes object it should be near
    #Calculate the colorbar position and size
    bbox = ax.get_position() 
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax  = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #Add the colorbar on the figure
    fig_object.colorbar(im, cax=cbar_ax)
    #Set label
    ax.set_ylabel(ylabel, labelpad = None, fontsize = 11)
    #disable axis
    ax.set_xticklabels([])




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



def TrimTrace(tr, trim_t1, trim_t2):
        #grab the signal and the noise
        #trace_trim    = tr.slice(starttime=start_window[0], endtime=end_window[1]).data
        trace_trim    = tr.slice(starttime=trim_t1, endtime=trim_t2)
        #Calculate RMS for signal and noise
        return trace_trim



#config_file = "configs.yaml"
#Open the configuration file
with open("confmtion.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Plot Option
Windowing      = Fig_params["Windowing"]
#Seimic mseed file
#mseedFiles        =   Fig_params['mseedFiles']
#Seimic mseed file
#Response_File     = Fig_params['Response_File']
#Grab the path of the data
#############################################
#Check the user want to make a plot bigger?
#plot_bigger    = Fig_params['plot_bigger']
# Plot the seismic Waveform 
waveform          = Fig_params['waveform']
pressure_waveform = Fig_params['pressure_waveform']
#Grab the bandpass
bandpass          = Fig_params['bandpass']
#get the start data 
STARTDATE          = Fig_params['STARTDATE']
########################################################
mseedFile_Land     = Fig_params["mseedFile_Land"]
######################################################
Response_File_Land = Fig_params["Response_File_Land"]
####################
#Land based station Response file

#get the window
start_window  = Fig_params['start_window']
end_window   = Fig_params['end_window']
#set the start and the endtime for signal
#startime_signal = Fig_params['startime_signal']
#endtime_signal  = Fig_params['endtime_signal']

#Set the number of subfigures to Zero
#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
#nfigs     = len(mseedFiles)
#Create the figure
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])


#periods List
pds     = bandpass.split("-")
f1      = 1./float(pds[-1]); f2 = 1./float(pds[-2]) ;
####################################################
f3      = 1./float(pds[1]);  f4 = 1./float(pds[0])
#f1 = 0.001; f2 = 0.00125 ; f3 = 0.002; f4 = 0.0025
pre_filt = [f1, f2, f3, f4]

#Create the figure
#Set the window
# Define signal and noise windows (replace with your own times)
#signal_window = (UTCDateTime("%sT%s"%(STARTDATE, start_window)), UTCDateTime("%sT%s"%(STARTDATE, end_window)))
t0_window = UTCDateTime("%sT%s"%(STARTDATE, start_window))
t1_window = UTCDateTime("%sT%s"%(STARTDATE, end_window))
#signal_window= (UTCDateTime("%sT%s"%(STARTDATE, startime_signal)), UTCDateTime("%sT%s"%(STARTDATE,endtime_signal)))


#TrimTrcae(tr, start_window, end_window)

fig, axs  = plt.subplots(3, 1, sharex = False, figsize = fig_size)
#fig, axs  = plt.subplots(1, 3, sharex = False, figsize = fig_size)
############## To check #################

#red the inventory
inv = read_inventory(Response_File_Land)
# Load waveform data
st = read(mseedFile_Land)  # Replace with your file path

#Filter using list comprehension, in order to keep only Z,N,E
#st = st.select(component=lambda x: x in ["Z", "N", "E"])
# Keep only channels ending in Z, N, or E
#st = st.select(channel="*Z") + st.select(channel="*N") + st.select(channel="*E")
#st = st.select(channel=lambda ch: ch[-1] in ["Z", "N", "E"])

##Remove instrument response
#st.remove_response(
#    inventory=inv,
#    output="VEL",         # Options: "DISP", "VEL", "ACC"
#    pre_filt=pre_filt,
#    water_level=60,       # Stabilizes inversion
#    zero_mean=True,
#    taper=True,
#    plot=False             # Optional: visualize frequency domain steps
#)


#Select components
#st_Z = st.select(component="Z")[0]
#st_N = st.select(component="N")[0]
#st_E = st.select(component="E")[0]
##################################

#parm_decov = "DISP"
parm_decov = "VEL"
st_Z = st.select(component="Z")[0].remove_response(inventory=inv, pre_filt=pre_filt, output=parm_decov,water_level=60)
st_N = st.select(component="N")[0].remove_response(inventory=inv, pre_filt=pre_filt, output=parm_decov,water_level=60)
st_E = st.select(component="E")[0].remove_response(inventory=inv, pre_filt=pre_filt, output=parm_decov,water_level=60)

##############################
st_Z =  st_Z.filter('lowpass', freq=1./50, corners=2, zerophase=True)
st_N =  st_N.filter('lowpass', freq=1./50, corners=2, zerophase=True)
st_E =  st_E.filter('lowpass', freq=1./50, corners=2, zerophase=True)

if(Windowing):
   st_Z = TrimTrace(st_Z, t0_window, t1_window)
   st_N = TrimTrace(st_N, t0_window, t1_window)
   st_E = TrimTrace(st_E, t0_window, t1_window)

##Define time window (adjust based on your data)
#start_time = st_Z.stats.starttime + 100
#end_time = start_time + 60
##Trim
#st_Z.trim(start_time, end_time)
#st_N.trim(start_time, end_time)
#st_E.trim(start_time, end_time)

# Filter (e.g., 0.1–1 Hz for surface waves)
#st_Z.filter("bandpass", freqmin=0.1, freqmax=1.0)
#st_N.filter("bandpass", freqmin=0.1, freqmax=1.0)
#st_E.filter("bandpass", freqmin=0.1, freqmax=1.0)

#st_Z.filter("bandpass", freqmin=0.0083, freqmax=0.1)
#st_N.filter("bandpass", freqmin=0.0083, freqmax=0.1)
#st_E.filter("bandpass", freqmin=0.0083, freqmax=0.1)
##################################################
#Get data arrays
z = st_Z.data * 1e+6
n = st_N.data * 1e+6
e = st_E.data * 1e+6
##Line with
lsize= 0.9
#Plot Z vs N (vertical–north)
axs[0].plot(n, z, lw = lsize, color='blue')
axs[0].set_xlabel("South-North (N)")
axs[0].set_ylabel("Z",  rotation=180, labelpad=15, fontsize=fsize)
#axs[0].set_ylabel("(Z)", fontsize=fsize)
#axs[0].set_title("Particle Motion (Z vs N)")
axs[0].grid(True)
#axs[0].axis("equal")
axs[0].set_ylim(-0.7, 0.5) #for  T>60s for the day 2024-03-29


#Optional: Plot Z vs E (vertical–east)
axs[1].plot(e, z, lw = lsize, color='green')
axs[1].set_xlabel("East-West (E)")
axs[1].set_ylabel("Z", rotation=180, labelpad=15, fontsize=fsize)
#axs[1].set_ylabel("(Z)", fontsize=fsize)
#axs[1].set_title("Particle Motion (Z vs E)")
axs[1].grid(True)
#axs[1].axis("equal")
axs[1].set_ylim(-0.7, 0.5) #for  T>60s for the day 2024-03-29

###############################
#Optional: Plot N vs E (vertical–east)
axs[2].plot(e, n, lw = lsize, color='k')
axs[2].set_xlabel("East-West (E)")
#axs[2].set_ylabel("(N)", fontsize=fsize)
axs[2].set_ylabel("N", rotation=180, labelpad=15, fontsize=fsize)
#axs[2].set_title("Particle Motion (N vs E)")
axs[2].grid(True)
#axs[2].axis("equal")
axs[2].set_ylim(-22, 21) #for  T>60s for the day 2024-03-29
#plt.show()
for i in range(len(axs)):
    axs[i].annotate(
        AphaDict[i],     # Text for the annotation
        #xy=(0.0089, 0.913),                   # Coordinates (relative to axes) for the annotation
        xy=xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=ftext,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
#axs[nfigs -1].set_xlabel('Time (hour:minute) on %s'%(DAYDATE), labelpad= pad, fontsize = fsize)
#figname   = "Waveforms_ParticleM.png"
if(Windowing):
    figname   = "Waveforms_ParticleM_windowing.pdf"
else:
    figname   = "Waveforms_ParticleM.pdf"
#figname   = "Waveforms_%s_of_%s.pdf"%(basename, DAYDATE)
#Set the font size of yticks
#plt.yticks(fontsize=fsize)
# Set the font size of xticks
#plt.xticks(fontsize=fsize)
#Make a rotation of xticks
#plt.xticks(rotation=45)
#Space between the subplots
#plt.subplots_adjust(hspace = 0.08)
plt.subplots_adjust(hspace = fig_space)
#Align ylabel of the figure
#fig.align_ylabels()
#Save the figure
for x_i in range(len(axs)):
    axs[x_i].tick_params(axis='both', labelsize=fsize)
    ###################
    #axs[x_i].tick_params(axis='y', labelrotation=90)
#figname = "TEST-PLOTS.png"
#fig.savefig(figname, bbox_inches = 'tight')
plt.savefig(figname, bbox_inches = 'tight')
#fig.savefig(figname)
