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


#xy_coords = (0.90, 0.85) # Coordinates (relative to axes) for the annotation
xy_coords = (0.0065, 0.84) # Coordinates (relative to axes) for the annotation

##Get the parameters
fsize = 14
pad   = 10
#Anotation fontsize
FsizeA=16
Figsize= 16
####################################

def ExtractorTemp(Filepath):
    df = pd.read_csv(Filepath)
    
    #Original columns
    ocols = [it for it in df.columns]
    
    #Arrange the data structure
    cols = [it.split('[')[0].strip() for it in df.columns]
    
    #Rename columns by replacing spaces in 'Battery Voltage' with a hyphen
    cols = ['%s-%s' % (it.split()[0], it.split()[1]) if 'Battery' in it else it for it in cols]
    
    try:
        df.rename(columns={ocols[1]: cols[1], ocols[2]: cols[2], ocols[3]: cols[3]}, inplace=True)
    except:
        raise ValueError("Not able to rename the columns %s  %s  %s" % (cols[1], cols[2], cols[3]))
    
    #Convert the time column into datetime64[ns]
    df['Time'] = pd.to_datetime(df['Time'], errors='coerce')
    
    #Transform the full temperature and time series into numpy arrays
    data_array = df['Temperature'].to_numpy()
    time_array = df['Time'].to_numpy()
    
    return (time_array, data_array)

####################################
#def ExtractorTemp(Filepath, DATEList):
#      df        =  pd.read_csv(Filepath)
#      #################################
#      #Originals columns
#      ocols     = [it for it in df.columns]
#      #Arrange the data structure
#      cols      = [it.split('[')[0].strip() for it in df.columns]
#      #change the colomn name by removig °, [ and switching Battery Voltage to Battery-Voltage
#      #Rename columns
#      cols      = ['%s-%s'%(it.split()[0], it.split()[1]) if 'Battery' in it else it  for it in cols]
#      try:
#           #RawData.data.rename(columns={ocols[1] : cols[1], ocols[2]: cols[2], ocols[3]: cols[3]}, inplace=True)
#           df.rename(columns={ocols[1] : cols[1], ocols[2]: cols[2], ocols[3]: cols[3]}, inplace=True)
#      except:
#            raise ValueError("Not able to renames the columns %s  %s  %s"%(cols[1], cols[2], cols[3])) 
#      #print(cols)
#      ##Convert the list of days to datetime.date objects
#      days_to_filter    = [pd.to_datetime(day).date() for day in list(DATEList)]
#      days_to_filter    = [pd.to_datetime(day).date() for day in days_to_filter]
#      #########################
#      #convert the time into  datetime64[ns]
#      df['Time']        = pd.to_datetime(df['Time'], errors='coerce')
#      #Filter rows where the date part of DateTime matches any of the days in the list
#      filtered_df       = df[df['Time'].dt.date.isin(days_to_filter)]
#      #Transform the Extracted Data into numpy array
#      data_array        = filtered_df['Temperature'].to_numpy()
#      time_array        = filtered_df['Time'].to_numpy()
#      #get the time and the data for sea level in m.asl
#      return(time_array, data_array)
#


def start_endtime(time_serie, ListDates):
    #minimum time
    tstart     = time_serie.min()
    #tstart     = time_serie[0]
    #print(tstart)
    #exit()
    #Add one day to the minimum time
    tend      = tstart + pd.to_timedelta(len(ListDates), unit= 'D')
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
def Plot_fig(ax, time, data, ListDates, param, labels, Plot_small= None, **PARAMSDICT):
    #get the parameters:
    #ylabel, color, ix = PARAMSDICT[param]
    ylabel, color = PARAMSDICT[param]
    #get the start and the endtime
    tstart, tend =  start_endtime(time, ListDates)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label = labels)
    if(Plot_small):
        #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.8, label = labels)
        ax.plot(time, data, lw=0.55, linestyle ="-", color = color, alpha =0.9, label = labels)
    else:
        ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label = labels)
        #ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0)
        #check if the user want to pot the Lake Level
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
    #if(float(ymax) < 0.1):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #get the coordinate of the point, where to write the text
#    #Write the text on the figure
#    ax.annotate(
#        AphaDict[ix],     # Text for the annotation
#        xy=xy_coords,                   # Coordinates (relative to axes) for the annotation
#        xycoords='axes fraction',    # Use axes fraction for positioning
#        fontsize=FsizeA,                 # Font size for the annotation
#        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
#        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
#    )
    #Add the legend
    ax.legend(loc="upper right")
    ##################
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #####################################
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #ax.legend(loc="upper right")
    #ax.legend(loc="upper left", fontsize='large',  bbox_to_anchor=(0, 1), frameon=True, shadow=False)
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #ax.yticks(fontsize = fsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
#    ax.set_ylim(-1.0, 1.0)
    #disable axis
    ax.set_xticklabels([])
    #Free memory by  manually trigger garbage collection
    #gc.collect()
##################********************************************************############################


#Define function
def PlotALL_fig(ax, time, data, param, labels, Plot_small= None, **PARAMSDICT):
    #get the parameters:
    #ylabel, color, ix = PARAMSDICT[param]
    ylabel, color = PARAMSDICT[param]
    #get the start and the endtime
    #tstart, tend =  start_endtime(time, ListDates)
    #Check the user need the plot to be bigger
    if(Plot_small):
        #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.8, label = labels)
        ax.plot(time, data, lw=0.35, linestyle ="-", color = color, alpha =0.8)
        #ax.plot(time, data, lw=0.55, linestyle ="-", color = color, alpha =0.9, label = labels)
    else:
        #ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label = labels)
        ax.plot(time, data, lw=1, linestyle ="-", color = color, alpha =0.8)
        #check if the user want to pot the Lake Level
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
    #if(float(ymax) < 0.1):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #get the coordinate of the point, where to write the text
#    #Write the text on the figure
    ax.annotate(
        labels,     # Text for the annotation
        xy=xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Add the legend
    #ax.legend(loc="upper right")
    ##################
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #####################################
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #ax.yticks(fontsize = fsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    #Set ylim of x-axis
#    ax.set_xlim(tstart, tend)
#    ax.set_ylim(-1.0, 1.0)
    #disable axis
    ax.set_xticklabels([])
    #Free memory by  manually trigger garbage collection
    #gc.collect()



















#config_file = "configs.yaml"
#Open the configuration file
with open("configALL.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Plot Option
Remove_frame    = Fig_params["Remove_frame"]
#Seimic mseed file
temperature     =   Fig_params['temperature']
#get the  temperature files
TEMPFILES       =   Fig_params['TEMPFILES']
#############################################
#DATELIST        =   Fig_params['DATELIST']
#Seimic mseed file
Response_File   = Fig_params['Response_File']
#Grab the path of the data
PARAMSDICT      = Fig_params['PARAMSDICT']
#############################################
#Check the user want to make a plot bigger?
plot_bigger    = Fig_params['plot_bigger']
###########
time_start    = Fig_params['time_start']
time_end      = Fig_params['time_end']
# Plot the seismic Waveform 
###########################
#Plot in a window, set this true if you're plottion a window in a one day
Windowing      = Fig_params['Windowing']
#Grab the bandpass
#Set the number of subfigures to Zero
#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs     = len(TEMPFILES)
#Create the figure
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

#ExtractorTemp(list(TEMPFILES)[0], list(DATELIST))
#exit()

#Create the figure
#fig, axs= plt.subplots(nfigs, 1, sharex = False, figsize=(12,10))

fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#####################################
#DATELIST  =  list(DATELIST)
#fixe the parameter
param     = "temperature"
#loop over the mseedFiles List
for ix,  TempFile in zip(range(len(TEMPFILES)), sorted(TEMPFILES)):
    #plot_bigger = False
    #print(TempFile)
    basename       = os.path.basename(TempFile)
    station        = basename.split(".")[0]
    #Extract date from the basename using regex
    #time, DataTemp = ExtractorTemp(TempFile, DATELIST)
    time, DataTemp = ExtractorTemp(TempFile)
    #print(DataTemp)
    #axs[ix].set_ylim([8.0, 12.0])
    axs[ix].set_ylim([7.0, 11.8])
    #exit()
    ###################################################
    #Plot the temperature
    #Plot_fig(axs[ix], time, DataTemp, DATELIST, param, station, **PARAMSDICT)
    PlotALL_fig(axs[ix], time, DataTemp, param, station, **PARAMSDICT)
    #print(DataTemp)

#get the DAYDATE
#DAYDATE = DATELIST[0]
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
        #ax.set_ylim([0.0, 10.5])
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
    #figname   = "Temp_%s_of_%s.pdf"%(basename, DAYDATE)
    figname   = "Temp_of_%s.pdf"%(basename)
    #format='%y %m %d %H %M %S'
    #axs[nfigs - 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[nfigs - 1].xaxis.set_major_formatter(mdates.DateFormatter("%m/%Y"))
#############################
#Plot a vertical line at a specific time.

#### Write on the x-axis
#axs[nfigs -1].set_xlabel('Time (hour:minute) on %s'%(DAYDATE), labelpad= pad, fontsize = fsize)
##########################################################
#Set the font size of xticks
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
