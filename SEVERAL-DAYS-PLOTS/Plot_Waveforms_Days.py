#!/usr/bin/env python

#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os,  xlrd, re, sys, math
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
import matplotlib.dates as mdates
import pandas as pd
from tabulate import tabulate
#for interpolation
#from scipy.interpolate import InterpolatedUnivariateSpline
#from datetime import datetime
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



AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            5: 'f', 6: 'f', 7:'h', 8:'g', 9: 'j',10:'k',
            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}



Figsize= 14
pad    =10
#############################
def check_if_string(var):
    return isinstance(var, str)
#################################
def inter1DD(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)




#Read the mseed file
def read_mseed(Date_of_Day, FILE, xmlFile, bandpass, Pressure):
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
#        data          = tr.data
        #time          = tr.times()
        #The sampling rate is 1 second
        date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= tr.data.size)
        #change the date_rng to numpy
        time      = date_rng.to_numpy()
        #time, data    = InterSpline(Date_of_Day, tr.times(), tr.data)
    #return(time, data, Title)
    return(time, tr, Title)





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

def start_endtime(time_serie, ListDates):
    #minimum time
    tstart     = time_serie.min()
    #tstart     = time_serie[0]
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




##Get the parameters
fsize = 14
pad   = 10

##############################################
#Define function
def Plot_fig(ax, time, data, ListDates, param, **PARAMSDICT):
    #get the parameters:
    ylabel, color, ix = PARAMSDICT[param]
    #get the start and the endtime
    tstart, tend =  start_endtime(time, ListDates)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label =param)
    ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label =param)
    ax.plot(time, data, lw=14, linestyle ="-", color = color, alpha =0.3, label =param)

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
        #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        xy=(0.0089, 0.87),                   # Coordinates (relative to axes) for the annotation
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
    #ax.legend(loc="upper right")
    #ax.legend(loc="upper left", fontsize='large',  bbox_to_anchor=(0, 1), frameon=True, shadow=False)
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #ax.yticks(fontsize = fsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
#    ax.set_ylim(-1.0, 1.0)
#    ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])

###############################################
def Extract_BS2D(df, Tlist,  param, nrange=None):
        if(nrange== None):
            r        = df.range.data
        else:
            r        = df.range.data[:nrange]
        if(param=="dlnEcho"):
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(df.amp, axis = 0)
                #Loop over the list of the time
                #mean_val = np.mean(np.mean(df_beam_avg, axis = 1), axis = 0)
                #Loop over the list of the time
                #Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist])
                Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist], axis=1)
                DList       = np.copy(Matrix2D)
                #for i in range(Matrix2D.shape[0]):
                for i in range(len(DList)):
                        #print(Matrix3D[beam_indx][i])
                        mean_val = np.mean(DList[i])
                        DList[i] = [(item -mean_val)/mean_val for item in DList[i]]
                #transform DList into an array
                #DList      = np.asarray(DList)[::-1]
                DList       = np.asarray(DList)
                #print(DList.shape)
                time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
                #return (time, r, Matrix3D[beam_indx])
                return (time, r, DList)



def Extract_ADCPData(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #get the range form the data
    r      = df.range.data
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
        T  = np.concatenate([df.velds.w.sel(time = ti)['time'] for ti in Tlist])
    elif(param=="Backscatter1D"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
    #    #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
    elif(param=="dlnEcho"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
        #Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist])
        Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist], axis= 1)
        T           = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
        DList = np.copy(Matrix2D)
        #for i in range(Matrix2D.shape[0]):
        #print(Matrix2D.shape)
        #exit()
        for i in range(len(DList)):
            mean_val = np.mean(DList[i])
            DList[i] = [(item -mean_val)/mean_val for item in DList[i]]
        DList = np.asarray(DList)
        ##############################
        P     = np.mean(DList, axis=0)

    else:
        #Loop over the  list of the time an extract the desire parameter
        P  = np.concatenate([df[param].sel(time = ti) for ti in Tlist])
        T  = np.concatenate([df[param].sel(time = ti)['time'] for ti in Tlist])
    #Return the  value corresponding to the date of your choice
    return (T, P)


def Plot_fig2D(ax, fig_object, df, time, ListDates, data2D, height_adcp, param, rframe= None, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #######################
    vmind = -40; vmaxd =+50
    ylabel, color, ix = PARAMSDICT[param] 
    #########################
    data2D = data2D * 100
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    #ylabel = "dln(Echo)";  
    ##color = "jet" 
    #color = "RdYlBu_r" 
    #############################
    ## Create custom colormap
    #ylabel = 'd'+r'$\ln$'+'(Echo)'  
    ###################################
    #vmind         =min(np.min(data2D, axis = 1))
    #vmaxd         = max(np.max(data2D, axis = 1))
    ########################
    #####Test sur L'ADCP orientation
    if(param=="ADCP_DOWN"):
        #invert the y-axis
        #invert ticks values from positive to negative
        depths    = abs((height_adcp+ blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
        y_lims    = [min(depths), height_adcp]
    else:
        depths    = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
        y_lims    = [height_adcp, max(depths)]
        #y_lims     = [max(depths), 0]
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, ListDates)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime'
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    #Set some generic y-limits.
    extent        = [start_num , end_num,  y_lims[0], y_lims[1]]
    #Make a 2D plot
    if(param=="ADCP_DOWN"):
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'upper', aspect = 'auto',
                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)

    else:
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'lower', aspect = 'auto',
                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax  = ax.get_ylim()
    xmin, xmax  = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #number of vertical lines for grid
    locations   = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks       = ax.get_yticks()
    ##################################
    if(param=="ADCP_DOWN"):
            #set the color bar on the figure
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.15, cbar_width, cbar_height])
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.1, 0.2, 0.3, 0.4])
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            cbar_ax.yaxis.set_label_position("right")
            #cbar_ax.ax.tick_params(labelsize=12)
            cbar_ax.yaxis.set_tick_params(labelsize=12)
            cbar_ax.yaxis.label.set_size(14)
            # Set spacing of 20 using MultipleLocator
            #cbar_ax.ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
            cbar_ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
            #position of the color bar
            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
            #ylabel and the postion of its positions
            ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="center", fontsize = Figsize)
            #ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
            ax.yaxis.set_label_coords(-0.08 , 1.2)
            #plt.yticks(fontsize = 12)
#        else:
#            #fig_object.colorbar(im, cax=cbar_ax)
#            ax.set_ylabel("H (m)", labelpad = pad, loc="top", fontsize = Figsize)
            #cbar_ax.yaxis.set_label_position("right")
    #Set label
    #Control font size
    ##############################################################
    # Set tick spacing to 10
    ax.yaxis.set_major_locator(MultipleLocator(5))
    #ax.set_yticklabels(ticks)
    ax.set_ylim(float(min(depths)), float(max(depths)))
    #Set anotations
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         xy=(0.0089, 0.93),                   # Coordinates (relative to axes) for the annotation
         #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=16,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    # Make ticks on x-axis and y-axis bigger
    #plt.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='major', labelsize=Figsize)
    #Remove the frame on the plot
    #ax.spines["right"].set_visible(False)
    if(rframe == "top"):
        ax.spines["top"].set_visible(False)
    elif(rframe == "bottom"):
        ax.spines["bottom"].set_visible(False)
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    elif(rframe == "right"):
        ax.spines["right"].set_visible(False)
    elif(rframe =="left"):
        ax.spines["left"].set_visible(False)
    #Set xlim of x-axis
#    ax.set_xlim(tstart, tend)
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])



def seismic_FFT(tr, Date_of_Day):
    #speed up figure plot
    #get the value from the Dictionary
    #get sampling rate
    sampling_rate = tr.stats.sampling_rate
    #get the data
    data          = tr.data
    #The sampling rate is 1 second
    #get the start and the endtime
    # Calculate spectrogram
    #freqs, times, spec = spectrogram(data, fs=sampling_rate, nperseg=256, noverlap=128, scaling='density')
    freqs, times, spec = scipy.signal.spectrogram(data, fs=1.0, nperseg=256, noverlap=128)
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






#config_file = "configs.yaml"
#Open the configuration file
with open("confsevdays.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Grab the Dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFiles     =   Fig_params['mseedFiles']
###########################################
#STARTDATE      =   Fig_params['STARTDATE']
#Seimic mseed file
Response_File       = Fig_params['Response_File']
#%%%%%%%%%%Pressure%%%%%%%%%%%%%%%%%%%%%%%%
Pressure            = Fig_params['Pressure']
#Grab the path of the data
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
##########################################
ADCP_DOWN           = Fig_params['ADCP_DOWN']
ADCP_UP             = Fig_params['ADCP_UP']
#ADCP parameter to plot
param_adpc          = Fig_params['param_adpc']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']
#############################################
#Check the user want to make a plot bigger?
#plot_bigger    = Fig_params['plot_bigger']
# Plot the seismic Waveform 
waveform       = Fig_params['waveform']
#######################
power_energy   =Fig_params['power_energy']
#########################
pressure_waveform= Fig_params['pressure_waveform']
#Grab the bandpass
bandpass       = Fig_params['bandpass']
#Set the number of subfigures to Zero
##change the values on the number of figure to plot by re-signing the len of PARAMSDICT
#nfigs     = len(mseedFiles)
##Create the figure
##Grab the space between the figure
#fig_space = Fig_params['fig_space']
##Grab the figure size
#fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])
##Create the figure
##fig, axs= plt.subplots(nfigs, 1, sharex = False, figsize=(12,10))
##fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#fig, ax  = plt.subplots(1, 1, sharex = False, figsize = fig_size)


#Set the number of subfigures to Zero
nfigs = 0
#Loop over the dictionary parameters
for i, key_param in zip(range(len(PARAMSDICT)), PARAMSDICT):
    #Check if the parameter is plotable
    if(Fig_params[key_param]):
        nfigs +=1
        #Add the index "i" in the list for numering the axis
        PARAMSDICT[key_param].append(i)

#Grab only the parameter that will be plotted by arranging the PARAMSDICT and re-asigned a new avalues to PARAMSDICT
PARAMSDICT = {k : PARAMSDICT[k]  for k in PARAMSDICT if(len(PARAMSDICT[k]) == 3)}
#Re-signed values to axis
for i, key in zip(range(len(PARAMSDICT)), PARAMSDICT):
        #try to get t he axis index
        PARAMSDICT[key][-1] = i

#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs  = len(PARAMSDICT)
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])
#Create the figure
#fig, axs= plt.subplots(nfigs, 1, sharex = False, figsize=(12,10))
fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)

####### Read the ADCP FILES##############
try:
    df_down = dlfn.read(ADCP_FILE_NAME_DOWN)
    df_up   = dlfn.read(ADCP_FILE_NAME_UP)
except:
    print("print provide the ADCP files:  %s and %s  "%(ADCP_FILE_NAME_UP, ADCP_FILE_NAME_DOWN))
#plot_discharge(ax1,hours_of_day, data_day, starttime, Remove_frame)

#Define a dictionary for Data storage
DictTraces = dict()
DictEnergy = dict()
#exit()
for ix,  mseedFile in zip(range(len(mseedFiles)), sorted(mseedFiles)):
    #plot_bigger = False
    basename = os.path.basename(mseedFile)
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
    time, tr, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure)
    #compute the seismic power
    tp, Power = seismic_FFT(tr, DAYDATE)
    #print(DAYDATE)
    if(waveform):
        DictTraces[DAYDATE]  = (time, tr)
    elif(power_energy):
        DictEnergy[DAYDATE] = (tp, Power)
    #Plot the figure by calling the plotting function, plot_twinx
    #time = time.astype("datetime64[ns]")
    #Plot_fig(axs[ix], ix, time, data, cmp_, waveform=True)
if(len(DictTraces) !=0):
    Edates                  = list(DictTraces.keys())
elif(len(DictEnergy) !=0):
    Edates                  = list(DictEnergy.keys())
#Extract the ADCP Data 
tad_up  ,   dlnBS_up    = Extract_ADCPData(df_up,   Edates,  param_adpc)      
tad_down,   dlnBS_down  = Extract_ADCPData(df_down, Edates,  param_adpc)      
#Convert pandas datetime to matplotlib date format
tadd       = mdates.date2num(tad_up)
##################################################
dlnBS_up   = dlnBS_up * 100
dlnBS_down = dlnBS_down * 100
#Extract the 2D Matrix perturbation for ADCP UP and DOWN
time_up,   rup, Matrix2D_up     = Extract_BS2D(df_up,  Edates,  param_adpc, nrange=None)
time_down, rdown, Matrix2D_down = Extract_BS2D(df_down,  Edates,  param_adpc, nrange=None)
#Check we need to plot the waveforms
if(waveform):
    Data  = np.concatenate([DictTraces[k][1].data for k in DictTraces.keys()])
    Time  = np.concatenate([DictTraces[k][0]      for k in DictTraces.keys()])
#Check if we need to plot seismic power
elif(power_energy):
    ####################
    #Concanete the Energy in the Dictionary
    eData  = np.concatenate([DictEnergy[k][1]      for k in DictEnergy.keys()])
    eTime  = np.concatenate([DictEnergy[k][0]      for k in DictEnergy.keys()])
    #interporlate the energy
    #tsnew, Powernew  =  inter1DD(ts,  Power, tadd)
    eTime, eData  =  inter1DD(eTime,  eData, tadd)
###################################
if(ADCP_UP):
    #set the parameter
    param      = "ADCP_UP"
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #plot the 2D Backscatter
    #plot the upward looking ADCP
    Plot_fig2D(axs[ix], fig, df_up, time_up, Edates, Matrix2D_up, ADCP_UP_HEIGTH_FROM_LAKEBED, param, rframe= None, **PARAMSDICT)

if(ADCP_DOWN):
    #set the parameter
    param      = "ADCP_DOWN"
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #plot the 2D Backscatter
    #plot the upward looking ADCP
    Plot_fig2D(axs[ix], fig, df_down, time_down, Edates, Matrix2D_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, param, rframe= None, **PARAMSDICT)
    
if(power_energy):
    param      = "power_energy"
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #Plot the seismic Energery
    Plot_fig(axs[ix], tad_up, eData, Edates, param, **PARAMSDICT)
    #Plot_fig(axs[ix], mdates.date2num(eTime), eData, Edates, param, **PARAMSDICT)
    #Create a second y-axis sharing the same x-axis
    ax2     = axs[ix].twinx()
    ##############
    #set a new color
    bs_color = 'r'
    ylabel_echo = 'd'+r'$\ln$'+'(Echo)'
    ###########################
    ax2.plot(tad_up, dlnBS_up, color= bs_color, lw=2, label='Echo')
    ax2.plot(tad_up, dlnBS_up, color= bs_color, lw= 14, alpha =0.3)
    #set the label
    ax2.set_ylabel(ylabel_echo, color=bs_color, labelpad = pad, fontsize=Figsize)
#print("****" * 10)
##Free memroy 
#Formatted the x-axis
if(len(Edates) ==1):
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
else:
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%B-%d"))
    #axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
plt.xticks(fontsize= 18)
##del DictTraces
#Save the figure
figname = "%s.png"%(basename)
fig.savefig(figname, bbox_inches = 'tight')
fig.savefig(figname)


##Set the axis-label for the last axis
##format='%y %m %d %H %M %S'
#axs[nfigs - 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
##axs[nfigs -1].xaxis.set_major_formatter(DateFormatter('%H:%M'))
#basename_ = os.path.basename(mseedFile) 
##### Write on the x-axis
#axs[nfigs -1].set_xlabel('Time (hour:minute) on %s'%(STARTDATE), labelpad= pad, fontsize = fsize)
##figname   = "Waveforms_%s_of_%s.png"%(basename, STARTDATE)
#figname   = "Waveforms_%s_of_%s.pdf"%(basename, STARTDATE)
##Set the font size of yticks
