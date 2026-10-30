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


#xy_coords = (0.0065, 0.87) # Coordinates (relative to axes) for the annotation
xy_coords = (0.0065, 0.86) # Coordinates (relative to axes) for the annotation

#Read the mseed file
def read_mseed(Date_of_Day, FILE, xmlFile, bandpass, Pressure = None, Lowpass=False):
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
        st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=30)
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
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
        #good for un
        if(Lowpass):
            tr            = tr.filter('lowpass', freq=0.1, corners=4, zerophase=True) #Good for Turbidity
            #tr            = tr.filter('lowpass', freq=1./25, corners=4, zerophase=True) #Good for Turbidity
        #sampling_rate = tr.stats.sampling_rate
        data          = tr.data
        #time          = tr.times()
        #The sampling rate is 1 second
        date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
        #change the date_rng to numpy
        time      = date_rng.to_numpy()
    return(time, data, Title)
    #return(time, tr, Title)


def ensure_list(parameter):
    if isinstance(parameter, list):
        return parameter
    return [parameter]



def Extract_database(df,  Dates, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    Tlist = ensure_list(Dates)
    df['date'] = df["DateTime"].dt.date
    dates_to_extract = [pd.to_datetime(year_month_day).date() for year_month_day in Tlist]
    #filtered for the selected days
    filtered = df[df['date'].isin(dates_to_extract)]
    #get the data for the selected days and the 
    # Step 4: Extract wind_speed
    Data_Select = filtered[['DateTime', param]]
    #get the time and the data values
    Time_ALL    = Data_Select['DateTime'].to_numpy()
    Data_ALL    = Data_Select[param].to_numpy()
    #Free the memory 
    del Data_Select, filtered
    #Return the  value corresponding to the date of your choice
    return (Time_ALL, Data_ALL)


def ExtractLakeLevel(df, Dates):
    Tlist = ensure_list(Dates)
    #Convert the list of days to datetime.date objects
    days_to_filter    = [pd.to_datetime(day).date() for day in Tlist]
    days_to_filter    = [pd.to_datetime(day).date() for day in days_to_filter]
    #convert the time into  datetime64[ns]
    df['DateTime']    = pd.to_datetime(df['DateTime'], errors='coerce')
    # Filter rows where the date part of DateTime matches any of the days in the list
    filtered_df       = df[df['DateTime'].dt.date.isin(days_to_filter)]
    #Transform the Extracted Data into numpy array
    data_array        = filtered_df['Data'].to_numpy()
    time_array        = filtered_df['DateTime'].to_numpy()
    #get the time and the data for sea level in m.asl
    return(time_array, data_array)



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

#if(waveform):
#ylabel = 'ACC (µm/s²)'
#ylabel = 'Pressure (Pa)'
##############################################
#Define function
def Plot_fig(ax, ix, time, data, ylabel, color, lw, Title=None):
    ########################################
    ylabel     = re.sub('', "", ylabel)
    #ax.plot(time, data, kwargs)
    #get the start and the endtime
    tstart, tend =  start_endtime(time)
    #Check the user need the plot to be bigger
    if(Title):
        network, station, component = Title.split("_")
        param      = "%s--%s"%(station, component)
        ##########
        ax.plot(time, data, lw=lw, linestyle ="-", color = color, alpha =0.9, label =param)
    else:
        ax.plot(time, data, lw=lw, linestyle ="-", color = color, alpha =0.9)

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
    # Add the legend ONLY if a label exists
    handles, labels = ax.get_legend_handles_labels()
    if labels:
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






#config_file = "configs.yaml"
#Open the configuration file
with open("confsev_meteo_lake.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Plot Option
Remove_frame      = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFiles        =   Fig_params['mseedFiles']
#STARTDATE       =   Fig_params['STARTDATE']
#Seimic mseed file
Response_File     = Fig_params['Response_File']
#Grab the path of the data
#############################################
#Check the user want to make a plot bigger?
#plot_bigger    = Fig_params['plot_bigger']
# Plot the seismic Waveform 
waveform          = Fig_params['waveform']
pressure_waveform = Fig_params['pressure_waveform']
###########################
#Plot in a window, set this true if you're plottion a window in a one day
Windowing         = Fig_params['Windowing']
#Grab the bandpass
bandpass          = Fig_params['bandpass']
#Set the number of subfigures to Zero
PlotDich    =  Fig_params['PlotDich']
PlotWind    =  Fig_params['PlotWind']
PlotSealevel= Fig_params['PlotSealevel']
#set a number to increase the axis
DictPlot = dict()
if(PlotWind):
    DictPlot[0]  = ['wind_speed', 'Wind speed (m/s)', 'teal']
    FILE_METEO = Fig_params['FILE_METEO']
    #read the Meteo database
    d_mteo = pd.read_csv(FILE_METEO, delimiter ='\t')
    #Add the of DateTime
    d_mteo = Add_Column(d_mteo)
if(PlotDich):
    FILE_DISCH = Fig_params['FILE_DISCH']
    #DictPlot['discharge']  = [1, 'Discharge (m³/s)','deepskyblue']
    if(PlotWind):
        DictPlot[1]  = ['discharge', 'Discharge (m³/s)','deepskyblue']
    else: 
        DictPlot[0]  = ['discharge', 'Discharge (m³/s)','deepskyblue']
    ###################################################
    d_disch                = pd.read_csv(FILE_DISCH, skiprows=0, sep =",",encoding="latin1")
if(PlotSealevel):
    FILE_LAKE = Fig_params['FILE_LAKE']
    if(PlotWind):
        if(PlotDich):
            DictPlot[2]  = ['Lake_Level', 'm a.s.l', 'teal']
    elif(PlotDich):
            DictPlot[1]  = ['Lake_Level', 'm a.s.l', 'teal']
    else:
            DictPlot[0]  = ['Lake_Level', 'm a.s.l', 'teal']
    #load the Lake Level Data
    dfLake = pd.read_csv(FILE_LAKE, skiprows=0, sep =",",encoding="latin1")


#print(d_disch) 

#Read the discharge data file
#Load the Discahge file
#Discharge file name
#print(DictPlot)
#Lake level File
#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs     = len(mseedFiles)
#Create the figure
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])


#Create the figure
if(len(DictPlot.keys())>0):
    nfigs_new = nfigs + len(DictPlot.keys())
    indexes = list(DictPlot.keys())
    start_index = len(indexes) 
    fig, axs  = plt.subplots(nfigs_new, 1, sharex = False, figsize = fig_size)
else:
    fig, axs    = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
    start_index = 0



#set the pressure variable here to be same nature a pressure_waveform
#Pressure     = pressure_waveform
color_signal = "k"
lw_signal =0.035
#loop over the mseedFiles List
#for ix,  mseedFile in zip(range(start_index, len(mseedFiles)), sorted(mseedFiles)):
#for ix, mseedFile in enumerate(sorted(mseedFiles), start=start_index):
for i, mseedFile in enumerate(sorted(mseedFiles)):
    #plot_bigger = False
    ix       =  start_index + i
    basename = os.path.basename(mseedFile)
    station = basename.split("_")[1]
    #Extract date from the basename using regex
    match    = re.search(r'\d{4}-\d{2}-\d{2}', basename)
    if(match):
        DAYDATE = match.group()
    else:
        print("No date for the file: %s"%(mseedFile))
        exit()
    ###################################################
    cmp_     = basename.split("-")[0]
    if("HH1" in cmp_):
        cmp_ = cmp_[:-1] +"Z"
    #get the seismogram from mseedFile
    if("HDH" in cmp_):
        time, data, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = True)
        ylabel = 'Pressure (Pa)'
        Plot_fig(axs[ix], ix, time, data, ylabel, color_signal, lw_signal, cmp_)
        #axs[ix].set_ylim(-5, 5) #for 2023-09-22
    else:
        #time, data, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = False)
        #time, data, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = False, Lowpass=True)
        time, data, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure = False, Lowpass=True)
        ####################
        data = data * 1e+6
        ###################
        ylabel_s = 'ACC (µm/s²)'
        Plot_fig(axs[ix], ix, time, data, ylabel_s, color_signal, lw_signal, cmp_)
        #if("HH1" in cmp_):
        #axs[ix].set_ylim(-1.0, 1.0) #Good for 2023-08-28
        #if("HHZ" in cmp_):
        #    axs[ix].set_ylim(-0.2, 0.2) #for 2023-09-22
        ##axs[ix].set_ylim(-5, 5) #for HH3 2023-09-22
#        axs[ix].set_ylim(-0.12, 0.12) #for HH3 2023-08-27
        ##Good for HH2
        #if(ix <=1):
        #    #axs[ix].set_ylim(-10, 10) #for 2023-09-22
        #    axs[ix].set_ylim(-1, 1)     #for 2023-08-27
        #else:
        #    #axs[ix].set_ylim(-5, 5) #for 2023-09-22
        #    axs[ix].set_ylim(-0.5, 0.5) #for 2023-08-27

##############################################
#Define function
#def Plot_fig(ax, ix, time, data, ylabel, color, lw, Title=None):

#get the seismic energy calculation
#########################
lw=2.5
if(PlotWind):
    idx    = indexes[0]  
    param  = DictPlot[idx][0]
    ylabel = DictPlot[idx][1]
    color  = DictPlot[idx][2]
    #read the Meteo database
    time, data         = Extract_database(d_mteo,  str(DAYDATE), param)
    #plot the data
    Plot_fig(axs[idx], idx, time, data, ylabel, color, lw)
if(PlotDich):
    ###################################################
    if(PlotWind):
        idx    = indexes[1]  
    else:
        idx    = indexes[0]  
    param  = DictPlot[idx][0]
    ylabel = DictPlot[idx][1]
    color  = DictPlot[idx][2]
    #read the Meteo database
    time, data         = ExtractLakeLevel(d_disch,  str(DAYDATE))
    #plot the data
    Plot_fig(axs[idx], idx, time, data, ylabel, color, lw)
if(PlotSealevel):
    #load the Lake Level Data
    if(PlotWind):
        if(PlotDich):
            idx    = indexes[2]
    elif(PlotDich):
        idx    = indexes[1]  
    else:
        idx    = indexes[0]  
    param  = DictPlot[idx][0]
    ylabel = DictPlot[idx][1]
    color  = DictPlot[idx][2]
    #read the Meteo database
    time, data         = ExtractLakeLevel(dfLake,  str(DAYDATE))
    #plot the data
    Plot_fig(axs[idx], idx, time, data, ylabel, color, lw)


#Set the font size of yticks
if(Windowing):
    ####################
    time_start   = Fig_params['time_start']
    time_end     = Fig_params['time_end']
    ######################################
    window_start = "%s %s"%(DAYDATE , time_start)
    window_end   = "%s %s"%(DAYDATE , time_end)
    time_1       = datetime.strptime(window_start, '%Y-%m-%d %H:%M')
    time_2       = datetime.strptime(window_end, '%Y-%m-%d %H:%M')
    #Create time series with 1-second sampling interval
    #print(time_1, time_2)
    #axs[-1].set_xlim([start, end])
    for ax in axs:
        #axs[ii].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
        #ax.axvline(specific_time_2, color='k', linestyle='--', alpha=0.4)
        #ax.set_xlim([start, end])
        ax.set_xlim([time_1, time_2])
#    #######################################################
#    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    window_start = window_start.replace(" ", "."); window_end = window_end.replace(" ", ".")
    figname      = "Waveform_%s_%s.pdf"%(window_start.replace(":", ".") ,window_end.replace(":", "."))
    #Rotate tick labels for better readability
    for label in axs[-1].get_xticklabels():
        label.set_rotation(45)
     #Set the axis-label for the last axis
else:
    basename_ = os.path.basename(mseedFile) 
    basename  = basename_.split(".mseed")[0]
    figname   = "Waveforms_%s_of_%s.pdf"%(basename, DAYDATE)
    #figname  = "Waveforms_%s_of_%s.png"%(basename, DAYDATE)
    #format='%y %m %d %H %M %S'
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#############################
#Plot a vertical line at a specific time.
#time_1 = pd.Timestamp('%s 09:27'%(DAYDATE))
#time_2 = pd.Timestamp('%s 11:52'%(DAYDATE))
##############################################
###time_2 = pd.Timestamp('%s 11:49'%(DAYDATE))
#time_3 = pd.Timestamp('%s 11:52:50'%(DAYDATE))
##time_2 = pd.Timestamp('%s 12:10'%(DAYDATE))
##axs[nfigs].axvline(specific_time, color='k', linestyle='--')
#for ax in axs:
#    #ax.axvline(time_1, color='r',lw=0.85, linestyle='--', alpha=0.3)
#    #ax.axvline(time_2, color='b',lw=0.85, linestyle='--', alpha=0.3)
#    ax.axvline(time_3, color='r',lw=0.85, linestyle='--', alpha=0.5)

##new cutt oof
########Good of  new zoom ########
#time_11 = pd.Timestamp('%s 07:50'%(DAYDATE))
time_11 = pd.Timestamp('%s 05:15'%(DAYDATE))
#time_11 = pd.Timestamp('%s 11:00'%(DAYDATE))
for ax in axs:
    ax.axvline(time_11, color='r',   lw=1.2, linestyle='--', alpha=0.8)
#    ax.axvline(time_22, color='r',lw=1.2, linestyle='--', alpha=0.8)
#    #ax.axvline(time_22, color='grey',lw=1.2, linestyle='--', alpha=0.95)
#

#### Write on the x-axis
axs[-1].set_xlabel('Time (hour:minute) on %s'%(DAYDATE), labelpad= pad, fontsize = fsize)
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
#fig.savefig(figname)
