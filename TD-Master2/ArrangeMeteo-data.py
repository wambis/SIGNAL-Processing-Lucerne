#!/usr/bin/env python

#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os,  xlrd, re, sys, math
import gc
import obspy
from obspy import read, read_inventory, Stream
from obspy.signal import PPSD
import matplotlib.pyplot as plt 
from glob import glob
from scipy import stats
from scipy.signal import detrend
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from math import log10, floor, ceil
import numpy as np
import matplotlib.dates as mdates
import pandas as pd
from tabulate import tabulate
#for interpolation
#from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import xarray as xr
import pandas as pd
from mhkit import dolfyn as dlfn
#################
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
#import matplotlib.dates as mdates
from mhkit.dolfyn.adp import api
from pandas.core.common import flatten
##################################
#from scipy.signal import spectrogram
import scipy
from scipy import interpolate
from obspy.core import UTCDateTime



#AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
#            5: 'e', 6: 'f', 7:'g', 8:'h', 9: 'j',10:'k',
#            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}



def check_if_string(var):
    return isinstance(var, str) 


def Extract_MeteoData(df,  Tlist, param):
    #Check if Tlist is a string
    if isinstance(Tlist, str):
        Tlist = [Tlist]
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
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



def Add_Column(ds):
    ds['DateTime'] = pd.to_datetime(ds[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
    #Return the dataframe
    return ds


def Extract_MeteoData(df,  Tlist, param):
    #Check if Tlist is a string
    if isinstance(Tlist, str):
        Tlist = [Tlist]
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def ExtractLakeLevel(df, Tlist):
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





#read the Meteo database
FILE_METEO = "DATA_ARRANGE_ALT23_24.csv"
#######################################


d_mteo     = pd.read_csv(FILE_METEO, delimiter ='\t')
#Add the of DateTime
d_mteo     = Add_Column(d_mteo)
print(d_mteo)
# Drop the unwanted columns
d_mteo.drop(columns=['Unnamed: 0', 'YR', 'MO', 'DA', 'HH', 'MM', 'SS'], inplace=True) 
# Move 'DateTime' to the first column
# Reorder columns with DateTime first
cols   = ['DateTime'] + [col for col in d_mteo.columns if col != 'DateTime']
d_mteo = d_mteo[cols]

# Optional: convert DateTime to datetime format if not already
d_mteo['DateTime'] = pd.to_datetime(d_mteo['DateTime'])

# Display the cleaned DataFrame
print(d_mteo.head())

# Save the cleaned DataFrame to a new file
d_mteo.to_csv('DataMeteo_cleaned_2023_2024.csv', index=False)
