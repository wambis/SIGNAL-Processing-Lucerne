#!/usr/bin/env python

#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os,  xlrd, re, sys, math
import gc
#import obspy
#from obspy import read, read_inventory, Stream
#from obspy.signal import PPSD
#from math import log10, floor, ceil
import matplotlib.pyplot as plt 
from glob import glob
from scipy import stats
from scipy.signal import detrend
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.dates as mdates
import pandas as pd
#for interpolation
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import pandas as pd
#################
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
##################################
#from scipy.signal import spectrogram
from obspy.core import UTCDateTime





def check_if_string(var):
    return isinstance(var, str) 


def  Extract_Data(df,  Tlist, param=None):
    #Check if Tlist is a string
    if isinstance(Tlist, str):
        Tlist         = [Tlist]
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Convert the list of days to datetime.date objects
    days_to_filter    = [pd.to_datetime(day).date() for day in Tlist]
    days_to_filter    = [pd.to_datetime(day).date() for day in days_to_filter]
    #convert the time into  datetime64[ns]
    df['DateTime']    = pd.to_datetime(df['DateTime'], errors='coerce')
    # Filter rows where the date part of DateTime matches any of the days in the list
    filtered_df       = df[df['DateTime'].dt.date.isin(days_to_filter)]
    #Transform the Extracted Data into numpy array
    if(param):
        data_array    = filtered_df[param].to_numpy()
    else:
        data_array    = filtered_df['Data'].to_numpy()
    #transform the time into an array
    time_array        = filtered_df['DateTime'].to_numpy()
    #get the time and the data 
    return(time_array, data_array)












#Grab the data for the river's discharge into the Lake Lucerne
FILE_DISCH  = "Discharge-Lucerne-2023-2024.csv"
#Load the Discahge file using pandas
d_disch     = pd.read_csv(FILE_DISCH, skiprows=0, sep =",",encoding="latin1")

#***************************************************************#
#Grab the meteorological data for a station located south of the Lake Lucerne
FILE_METEO = "DataMeteo_cleaned_2023_2024.csv"
#Load the meteorological file using pandas
d_mteo     = pd.read_csv(FILE_METEO, skiprows=0, sep =",",encoding="latin1")

#***************************************************************#
#Grab the data for the level of the Lake Lucerne
FILE_LAKE  = "Lake-Level-Lucerne-2023-2024.csv"
#Load the  Lake Level file using pandas
dfLake      = pd.read_csv(FILE_LAKE, skiprows=0, sep =",",encoding="latin1")


#Plot the Data for this selected date
SELECDATE = ["2023-10-23"]





#Define the figure
fig, axs = plt.subplots(3, 1, figsize=(12, 10))


#Extract the wind speed on the 2023-10-23
#set the parameter to plot, "wind_speed" as defined in the data file
param                   = "wind_speed"
#Call the the function ``Extract_Data`` to extract the desire data
time_wind, data_wind    = Extract_Data(d_mteo,  SELECDATE, param)


#Extract the Lake Level on the 2023-10-23
#Call the the function ``Extract_Data`` to extract the desire data
time_level, data_level = Extract_Data(dfLake,  SELECDATE)


#Extract the River Dischage on the 2023-10-23

#Call the the function ``Extract_Data`` to extract the desire data
time_d, data_d     = Extract_Data(d_disch,  SELECDATE)

#Plot the data on the figure, note that the figure has 3 axes, namely axs[0], axs[1], axs[2]
#plot the wind speed on the first axis, axs[0]
axs[0].plot(time_wind, data_wind, label='wind speed', color='k', lw=3)
axs[0].legend()
axs[0].grid(True)
#disable x-axis ticks
axs[0].set_xticklabels([])
#Set ylabel
axs[0].set_ylabel('wind speed (m/s)' , fontsize = 12)

#plot the Lake Level on the second axis, axs[1]
axs[1].plot(time_level, data_level, label='lake level', color='blue', lw=3)
axs[1].legend()
axs[1].grid(True)
#Set y-axis to meters above sea level
axs[1].get_yaxis().get_major_formatter().set_useOffset(False)
axs[1].get_yaxis().set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.2f}'.format(val)))
#disable x-axis ticks
axs[1].set_xticklabels([])
#Set ylabel
axs[1].set_ylabel('lake level (m.a.s.l)' , fontsize = 12)

#plot the river discharge on the third axis, axs[2]
axs[2].plot(time_d, data_d, label='discharge', color='lightblue', lw=3)
axs[2].legend()
axs[2].grid(True)
#Set ylabel
axs[2].set_ylabel('Discharge (m³/s)' , fontsize = 12)

#Format the last axis
axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
#set the label of the last axis
#label the x-label
#axs[-1].set_xlabel('Time (hour:minute) on %s'%(SELECDATE[0]), labelpad =12, fontsize = 11)
axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

plt.xticks(fontsize= 18)
########################################
#make sure all the subplots are aligned
fig.align_ylabels(axs)
#set the figure name
figname = "Figure-Example.png"
#Save the figure
fig.savefig(figname, bbox_inches = 'tight')
