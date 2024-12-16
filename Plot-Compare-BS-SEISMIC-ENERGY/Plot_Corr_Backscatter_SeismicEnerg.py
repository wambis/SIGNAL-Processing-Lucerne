#!/usr/bin/env python

#Plot the mseed seismic data files
import obspy
from obspy import read, read_inventory
import matplotlib.pyplot as plt 
from glob import glob
import numpy as np
#from Ipython.display import display
import pandas as pd
#for interpolation
import xlrd
from obspy import read, read_inventory
from glob import glob
from scipy import stats
import matplotlib.ticker as mticker
from math import log10, floor, ceil
from tabulate import tabulate
#for interpolation
from scipy import interpolate
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import matplotlib.dates as dt
import xarray as xr
from mhkit import dolfyn as dlfn
from mhkit.dolfyn.adp import api
from datetime import datetime
#################
from scipy.fft import fft, fftfreq
#from scipy.signal import spectrogram
import scipy
##############
import matplotlib.dates as mdates
##################################
from obspy.core import UTCDateTime
#############################################
from sklearn.linear_model import LinearRegression 
from sklearn.metrics import mean_squared_error, r2_score 
#####################################################
from sklearn.model_selection import train_test_split 
##################################################
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, MinMaxScaler, StandardScaler 
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
import datetime


def inter1D(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

def inter1BS(x, y):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    t24h =  24 * 60
    xnew = np.arange(0, t24h, 10)
    #exit()
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

def inter1DD(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    #t24h =  24 * 60
    #xnew = np.arange(0, t24h, 10)
    #start = 0 
    #stop = 3600 * 24 
    #num_values = 144
    #step = (stop - start) / (num_values - 1)
    #xnew = np.arange(min(x), stop + step, step)
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    #print(x)
    #print("****" * 20)
    #print(xnew)
    #exit()
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

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


#Function to create and fit the polynomial model
def fit_polynomial(degree, X, y):
    #Create a Polynomial pipeline
    model = make_pipeline(PolynomialFeatures(degree), LinearRegression())
    model.fit(X,y)
    return model

#def read_mseed(FILE):
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
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",water_level=60)
    #get the trace from the stram
    for tr in st:
        #get the parameters of the trace
        network       = tr.stats.network
        station       = tr.stats.station
        channel       = tr.stats.channel
        stime         = str(tr.stats.starttime)
        endtime       = str(tr.stats.endtime)
        #remove mean and trend on the trace
        tr.detrend(type='demean')
        #Title         = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
        Title         = "%s  %s"%(network, station)
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
    return(time, tr, Title)


def seismic_FFT(tr, Date_of_Day):
    #speed up figure plot
    #get the value from the Dictionary
    #get sampling rate
    sampling_rate = tr.stats.sampling_rate
    #get the data
    data          = tr.data
    #The sampling rate is 1 second
#    date_rng      = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
    #change the date_rng to numpy
#    time          = date_rng.to_numpy()
    #get the start and the endtime
    # Calculate spectrogram
    #freqs, times, spec = spectrogram(data, fs=sampling_rate, nperseg=256, noverlap=128, scaling='density')
    freqs, times, spec = scipy.signal.spectrogram(data, fs=1.0, nperseg=256, noverlap=128)
    #freqs, times, spec = scipy.signal.spectrogram(data, fs=600, nperseg=256, noverlap=128)
    ##########
    #which axes object it should be near Calculate the colorbar position and size
    #set the colorbar annotation
    # Calculate PSD 
    #################################
    #ENER_           = np.mean(spec, axis =0)
    ENER_           = 10 * np.log10(np.mean(spec, axis =0))
    #ENER_           = 10 * np.log10(np.mean(spec, axis =1))
    # Calculate PSD PSD = np.abs(Sxx)**
    #Call the interpolation
#    tf, tc           = inter1D(times, times, freqs)
#    #time_new, freqs_new = inter1D(times[:len(freqs)], freqs, time)
#    time_new, freqs_new = inter1D(tf, freqs, time)
#    return(time, freqs_new)
    return(times, ENER_)


def plot_spectrogram(ax,tr, Date_of_Day):
    #speed up figure plot
    #get sampling rate
    sampling_rate = tr.stats.sampling_rate
    #get the data
    data          = tr.data
    #The sampling rate is 1 second
    date_rng      = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
    #change the date_rng to numpy
    time          = date_rng.to_numpy()
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime' 
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    ##############################33
    extent        = [start_num , end_num]
    #Plot the spectrogram
    spec, freqs, times, img = ax.specgram(data, NFFT=256, noverlap = 128, Fs=sampling_rate,
                                        #detrend='mean', mode ='psd', cmap = 'RdYlBu', xextent = extent)
                                        detrend='mean',  cmap = 'RdYlBu_r', xextent = extent)

    ##########
    ## Define the boundaries for the discrete intervals
    #levels = MaxNLocator(nbins=4).tick_values(spec.min(), spec.max())
    ## Create a BoundaryNorm object
    #norm = BoundaryNorm(levels, ncolors=256, clip=True)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    image_data = img.get_array()
    PSD_MEANS_t  = np.mean(image_data, axis = 0)
    PSD_MEANS_freq  = np.mean(image_data, axis = 1)
    #################################
    #ENER_     = np.sum(np.square(spec))
    ENER_      = np.mean(spec, axis =0)
    #ENER_      = np.square(spec)
    ENRER_db   = 10 * np.log10(ENER_)
    #print(spec.shape)
    #print(len(freqs), len(times))
    #print("***"*10) 
    #print(image_data.shape) 
    #print(len(PSD_MEANS_freq), len(PSD_MEANS_t))
#    print("-----"*50)
#    print(PSD_MEANS)
#    print("-----"*50)
#    print(len(ENRER_db), len(PSD_MEANS))
    #file_name  = 'SEMIC_PDS_%s.dat'%(Date_of_Day)
    file_name  = 'SEMIC_ENERGY_%s.dat'%(Date_of_Day)
    #Save the array into na file
    np.savetxt(file_name , PSD_MEANS_t , delimiter="\n", fmt="%.4f")
    #np.savetxt(file_name , ENRER_db , delimiter="\n", fmt="%.4f")



def check_if_string(var): 
    return isinstance(var, str)
#Defined function for data extraction

def Extract_ADCPData(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
    if(check_if_string(Tlist)):
        Tlist =[ Tlist ]
    #Check if the desired velocity is the Horizontal velocity
    if(param=="velocity"):
        #Loop over the list of the time
        try:
            P  = np.concatenate([np.nanmean(df.velds.U_mag.sel(time = ti), axis =0) for ti in Tlist])
            T  = np.concatenate([df.velds.U_mag.sel(time = ti)['time'] for ti in Tlist])
        except:
            print("The Data seems not to exist for the certain data in the list Time :  ")
            print(' , '.join(Tlist))
            #exit()
            #P  = np.concatenate([np.nanmean(df.velds.U_mag.sel(time = ti, method='nearest'), axis =0) for ti in Tlist])
            #T  = np.concatenate([df.velds.U_mag.sel(time = ti,method='nearest')['time'] for ti in Tlist])

    #Check if the desired velocity is the vertical velocity
    elif(param=="vertical_vel"):
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df.velds.w.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df.velds.w.sel(time = ti)['time'] for ti in Tlist])
    elif(param=="Backscatter1D"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0) 
    #    #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])

    else:
        #Loop over the  list of the time an extract the desire parameter
        P  = np.concatenate([df[param].sel(time = ti) for ti in Tlist])
        T  = np.concatenate([df[param].sel(time = ti)['time'] for ti in Tlist])
    #Return the  value corresponding to the date of your choice
    return (T, P)


def Add_Column(ds):
    ds['DateTime'] = pd.to_datetime(ds[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
    #Return the dataframe
    return ds

#
def Average(df, n_bin=1):
    #To average the data into time bins, we need the parmeter
    #Define functio
    #1) n_bin which represents the number of data points in each ensemble,
    #2) here our enseblme happen every 600s and  the frequency is fs = 0.001668, so n_bin = fs* 600 = 1
    #n_bin = 1 is the value by default
    #Set the average tool
    avg_tool = dlfn.VelBinner(n_bin = n_bin, fs = df.fs)
    #Perform the averge
    df_avg   = avg_tool.bin_average(df)
    #return df_avg
    return df_avg



def is_sorted(arr):
    for i in range(len(arr) - 1):
        if arr[i] > arr[i + 1]:
            return False
    return True

## Example usage
#arr = [1, 2, 3, 4, 5]
#print(is_sorted(arr))  # Output: True



#set the date
Date_of_Day = "2023-10-23" 
#Date_of_Day = "2023-10-20" 
########*********************##################
xmlFile     = 'RESPOSE-FILE.XML/MUS03_response.xml'
mseedFile   = 'DATASELECT/STATION03/XJ_MUS03_HH2-2023-10-23.mseed'
#mseedFile   = 'DATASELECT/STATION03/XJ_MUS03_HH2-2023-10-20.mseed'
#mseedFile   = 'DATASELECT/STATION03/XJ_MUS03_HH3-2023-10-23.mseed'
bandpass    = '0.02-0.022-200-1000'
######################################
ADCP_FILE_NAME = 'DATASELECT/MUO2_DOWN_11_23.000'
#ADCP_FILE_NAME = 'DATASELECT/MUO2_UP_11_23.000'

#number bin
n_bin      = 1
#***Just for pressure**********#
df_adcp    = dlfn.read(ADCP_FILE_NAME)
df_avg     = Average(df_adcp, n_bin)
#beam_up  = df_up_av.beam_angle
#blank_up = df_up_av.blank_dist
#####*****************#######
#get the beam angle
beam_angle   = df_avg.beam_angle
#get the blink distance
blank_dist   = df_avg.blank_dist
############################################
ranges       = df_avg.range.data
#set the parameter
param        = "Backscatter1D"
#Extract the Data
tad  , d_amp = Extract_ADCPData(df_avg,  Date_of_Day, param)
# Convert pandas datetime to matplotlib date format 
tadd    = mdates.date2num(tad)
#####get the seismogram from mseedFile
time, tr, Title   = read_mseed(Date_of_Day, mseedFile, xmlFile, bandpass)
##
ts, Power         = seismic_FFT(tr, Date_of_Day)
##print(tad) 
#Power             = np.loadtxt("SEMIC_ENERGY_2023-10-23.dat") 
#Power             = np.loadtxt("SEMIC_ENERGY_2023-10-23_meantime.dat") 
#Power =  np.loadtxt("SEMIC_ENERGY_2023-10-23_meansFreqs.dat") 

#d_amp = np.asarray(sorted(d_amp))
##############################
#Power    = np.asarray(sorted(Power))
#tsnew, Powernew =  inter1D(ts,  Power, tadd)

#Perform the interpolation
tbsnew, d_ampnew  =  inter1D(tadd,  d_amp, Power)
#
#wamba
#We sorted the backscatter because it is a proxy of sediment load and consider that the 
#sediment crease progressively in the lake, that the means the quantity of the sediment increase
#as function of time
d_ampnew      = np.array(sorted(d_ampnew))
#d_amp      = np.array(sorted(d_amp))
#set parameters
predict = True
pad = 10
fsize = 14
polynomial_degree = 10 
#polynomial_degree = 9 
#polynomial_degree = 7 
#polynomial_degree = 6 
#
###set the figure size
fig_width  = 12
fig_height = 11
#
#Grab the figure size
fig_size   = (fig_width, fig_height)
fig, ax    = plt.subplots(1, 1, sharex = False, figsize = fig_size)
#plot the frequency as function of time

#set value for training
X = d_ampnew.reshape(-1, 1)
y = Power

#######################################

#Split Data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
#Polynomial Features

####Pipeline
pipeline  = Pipeline([('poly_features', PolynomialFeatures(degree=polynomial_degree, include_bias=False)), 
                     ('scaler', StandardScaler()), 
                     ('linear_regression', LinearRegression())])

#Train the model 
pipeline.fit(X_train, y_train)

#Predict on the testing data 
y_pred = pipeline.predict(X_test)

#Evaluate the model 
mse = mean_squared_error(y_test, y_pred) 
print(f'Mean Squared Error: {mse}')

#
### Generate polynomial features 
## Visualize the results 
if(predict):
    #plot the original data
    ax.scatter(X, y, color='k', label='Original data', s=22, edgecolors='white', alpha = 0.98) 
    #make the x-values to plot the predicted values
    X_plot = np.linspace(X.min(), X.max(), X.shape[0]).reshape(-1, 1) 
    y_plot = pipeline.predict(X_plot) 
    #plot the predicted values
    ax.scatter(X_plot, y_plot, label='Polynomial Fit', color='red', s= 8, edgecolors='m', alpha= 0.7)
    #Add the legend
    ax.legend(loc="upper right",fontsize=fsize )
    figname = "Fig_Coor_PLOT_%s_predict.png"%(Date_of_Day)
else:
    #ax.scatter(d_amp, y, color='k',s=10 )
    ax.scatter(X, y, color='k', label='Original data',s=22, edgecolors='white', alpha = 0.98) 
    #Add the legend
    ax.legend(loc="upper right",fontsize=fsize )
    figname = "Fig_Coor_PLOT_%s.png"%(Date_of_Day)
##################
#Add the legend
ax.set_title('Seismic energy as a function of backscatter on %s' %(Date_of_Day), fontsize=fsize, loc ='center', pad =pad)
#Make the grid of the axis
ax.grid(visible = True, axis = "both", alpha = 0.9)
#Set label
ax.set_ylabel("Seismic energy (dB)", labelpad = pad, fontsize = fsize)
ax.set_xlabel("Backscatter (counts)", labelpad = pad, fontsize = fsize)
plt.yticks(fontsize = 11)
# Make ticks on x-axis and y-axis bigger 
plt.tick_params(axis='both', which='major', labelsize=14) 
#Save figure 
plt.savefig(figname, bbox_inches = 'tight')
