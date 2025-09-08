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


xy_coords = (0.0065, 0.89) # Coordinates (relative to axes) for the annotation

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
###########################################################################

from obspy import read
from obspy.signal.polarization import particle_motion_odr
import matplotlib.pyplot as plt

# Load three-component waveform data
st = read("your_data.mseed")  # Replace with your file path

# Sort the stream to ZNE order (important!)
st.sort(['channel'])

# Trim to Rayleigh wave window
start = st[0].stats.starttime + 100
end = start + 60
st.trim(start, end)

# Apply bandpass filter to isolate surface waves
st.filter("bandpass", freqmin=0.1, freqmax=1.0)

# Run particle motion analysis
azimuth, incidence, azimuth_err, incidence_err = particle_motion_odr(st)

# Print results
print(f"Azimuth: {azimuth:.2f}° ± {azimuth_err:.2f}")
print(f"Incidence: {incidence:.2f}° ± {incidence_err:.2f}")

#############################################################
##############################################################
z = st.select(component="Z")[0].data
n = st.select(component="N")[0].data
e = st.select(component="E")[0].data
###########################################
# Plot Z vs N
plt.figure(figsize=(6, 6))
plt.plot(n, z, color='blue')
plt.xlabel("North-South (N)")
plt.ylabel("Vertical (Z)")
plt.title("Rayleigh Wave Particle Motion (Z vs N)")
plt.grid(True)
plt.axis("equal")
plt.show()




