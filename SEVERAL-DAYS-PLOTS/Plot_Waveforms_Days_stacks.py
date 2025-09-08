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


AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            5: 'f', 6: 'g', 7:'g', 8:'i', 9: 'j',10:'k',
            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}

##set the coordinate
#xy_coords = (0.0065, 0.953) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.93) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.91) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.90) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.89) # Coordinates (relative to axes) for the annotation
xy_coords = (0.0065, 0.85) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.82) # Coordinates (relative to axes) for the annotation
#xy_coords = (0.0065, 0.73) # Coordinates (relative to axes) for the annotation

#Figsize= 14
#Figsize= 18
Figsize= 16
pad    = 20
#Anotation fontsize
FsizeA=20
#FsizeA=18
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

def inter1D(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

#Define the bins calculation 
def hist_nbin(data):
    nbins = int(1+3.22 * np.log10(data.size))
    #nbins = np.arange(0, nbins +1, 0.15)
    nbins = np.arange(0, nbins +1, 1)
    return(nbins)


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


##############################################
#Define function
def Plot_fig(ax, time, data, ListDates, param, Plot_small= None, **PARAMSDICT):
    #get the parameters:
    ylabel, color, ix = PARAMSDICT[param]
    #get the start and the endtime
    tstart, tend =  start_endtime(time, ListDates)
    #Check the user need the plot to be bigger
    #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
    #ax.plot(time, data/max(data), lw=0.45, linestyle ="-", color = color, alpha =0.9, label =param)
    if(Plot_small):
        ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.8, label =param)
        #ax.plot(time, data, lw=0.8, linestyle ="-", color = color, alpha =0.9, label =param)
    else:
        ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label =param)
        ax.plot(time, data, lw=12, linestyle ="-", color = color, alpha =0.5, label =param)
        #check if the user want to pot the Lake Level
    #if(param == "Lake_Level"):
    if(param == "Lake_Level" or param=="pressure_up"):
        #Scatter the data first
        #Set y-axis to meters above sea level
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.2f}'.format(val)))
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
    #if(float(ymax) < 0.1):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #if (gridloc_rmv==None):
    #    #number of vertical lines for grid
    #    locations = verical_grid()
    #    ax.xaxis.set_major_locator(locations)
    #    #Make the grid of the axis
    #    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #else:
    #    print("No grid is chooseen")
    #Add the legend
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
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
#    ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])
    #Free memory by  manually trigger garbage collection
    #gc.collect()

###plot the envelope
def PlotEnvelop(ax, time,  tr_filt, ListDates, param, Pressure= False, linewidth=1.2, **PARAMSDICT):
    #get the sampling rate
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    #sampling_rate = tr_filt.stats.sampling_rate
    #check if the data is pressure, then change the ylabel to Pa
    if(Pressure):
        ylabel = "Pa"
    #get the data from the filtered trace
    #tr_filt.filter('bandpass', freqmin=1.0/200, freqmax=1.0/167, corners=2, zerophase=True) #Good for 2023-10-23
    tr_filt.filter('bandpass', freqmin=1.0/200, freqmax= 1.0/158, corners=2, zerophase=True) #Good for 2023-10-23
    #tr_filt.filter('bandpass', freqmin=0.006, freqmax= 0.015, corners=2, zerophase=True) #Good for 2023-09-02
#    tr_filt.filter('bandpass', freqmin=0.004, freqmax= 0.01, corners=2, zerophase=True) #Good for 2023-09-02
#    tr_filt.filter('bandpass', freqmin=0.002, freqmax= 0.009, corners=2, zerophase=True) #Good for 2023-09-02
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, ListDates)
    #Envelope of filtered data
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    #Plot the amplitude spectrum
    #ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5, label = AphaDict[ix])
    ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5)
    #Fill the space between the Figure and the x-xais
    #ax.fill_between(time, data_envelope, color = 'r', alpha = 0.3)
    ax.fill_between(time, data_envelope, color = color, alpha = 0.3)
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        #xy=(0.0089, 0.94),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.92),                   # Coordinates (relative to axes) for the annotation
        xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize= FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Remove the frame on the plot
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ##################################
    #if(float(ymax) < 0.01):
#        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    #ax.legend(loc="upper left", fontsize='large', bbox_to_anchor=(0, 1), frameon=True, shadow=False)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
#    #Make the grid of the axis
#    ax.grid(visible = True, axis = "x", alpha = 0.7)
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
    #Set label
    #disable ticks on  x-axis
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    #Remove the ticks marks (small vertical lines) on the x-axis
    ax.tick_params(axis="both", length = 0, color= "white", width = 0)
    ax.tick_params(axis='both', labelsize=Figsize)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #Free memory by  manually trigger garbage collection
    gc.collect()

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






def control_cbar_label(cbar):
    #Control font size
    cbar.yaxis.set_tick_params(labelsize=12)
    cbar.yaxis.label.set_size(14)
    #Set spacing of 20 using MultipleLocator
    #cbar_ax.ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
    cbar.yaxis.set_major_locator(mticker.MultipleLocator(20))
    cbar.tick_params(labelsize=10)


#%%%%%%%%%%%%%%%%%%%%%%%%
def Plot_fig2D(ax, fig_object, df, time, data2D, ListDates, height_adcp,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    ylabel, color, ix  = PARAMSDICT[param]
    vmind         = min(np.min(data2D, axis = 1))
    vmaxd         = max(np.max(data2D, axis = 1))
    #if("velocity" in ylabel):
    ##if("velocity2D_" in param):
    #    vmind = 0.0; vmaxd = 0.04
    if("dlnEcho" in param):
        #vmind = -50; vmaxd =+50
        vmind = -40; vmaxd =+40
    #####Test sur L'ADCP orientation
    if(ADCP_DOWN):
        #invert ticks values from positive to negative
        depths     = abs((height_adcp+ blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
        y_lims     = [min(depths), height_adcp]
    else:
        depths     = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
        y_lims     = [height_adcp, max(depths)]
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
    if(ADCP_DOWN):
        im= ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'upper', aspect = 'auto',
        interpolation ='bilinear', interpolation_stage='rgba', resample=True, extent = extent)
        #interpolation ='spline16', resample=True, extent = extent)
        #set the limit
        ax.set_ylim(min(depths), float(height_adcp))
    else:
        im = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'lower', aspect = 'auto',
        interpolation ='bilinear',interpolation_stage='rgba', resample=True, extent = extent)
        #interpolation ='spline16', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    #cbar_x      = bbox.x1 + 0.09
    cbar_y      = bbox.y0
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    ##################################
    ##################################
    if('CW' in ylabel):
        #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
        if(ADCP_DOWN):
            #set the color bar on the figure
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.05, cbar_width, cbar_height])
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
            cbar_ax.yaxis.set_label_position("right")
            #shift the label by 3.5 along x-axis and by 1.0 along the y-axis
            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
            #move the label by -0.08 on x-axis and by 1.3 on y-axis
            #ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = Figsize)
            #ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            ax.set_ylabel("H (m)", labelpad = pad, fontsize = Figsize)
            ##################################################
            #Change the position of the y-label, here we juste changed the y-poistion, while keeping the x-position fix
            label_y = ax.yaxis.get_label()
            x, _  = label_y.get_position()
            label_y.set_position((x, 1.3))
            #Control font size
            cbar_ax.yaxis.set_tick_params(labelsize=12)
            cbar_ax.yaxis.label.set_size(14)

            #Set font size for the label and ticks
            cbar_ax.tick_params(labelsize=Figsize)
            # Set spacing of 20 using MultipleLocator
            #control the color bar label
            #control the color bar label
            #control_cbar_label(cbar_ax)
    else:
        if(ADCP_DOWN):
            #set the color bar on the figure
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.05, cbar_width, cbar_height])
            #cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.07, cbar_width, cbar_height])
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.1, 0.2, 0.3, 0.4])
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.01, 0.02, 0.03, 0.04])
            cbar_ax.yaxis.set_label_position("right")
            #ylabel and the postion of its positions
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
#            x_position = -0.07 - pad / 1000  # Adjust scaling factor as needed
#            ax.yaxis.set_label_coords(x_position , 1.3)
            ax.set_ylabel("H (m)", labelpad = pad , loc="center", fontsize = Figsize)
            #Change the position of the y-label, here we juste changed the y-poistion, while keeping the x-position fix
            label_y = ax.yaxis.get_label()
            x, _  = label_y.get_position()
            label_y.set_position((x, 1.3))
            #Set font size for the label and ticks
            cbar_ax.tick_params(labelsize=Figsize)
            #control the color bar label
            #Set font size for the label and ticks
            ##Now set font size and padding through the colorbar's axis
            #cbar_ax.yaxis.label.set_size(Figsize)        # Label font size
            #cbar_ax.yaxis.labelpad = pad                  # Padding
#######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(not ADCP_DOWN):
        #ax.set_ylim(float(min(depths)), float(max(depths)))
        ax.set_ylim(float(height_adcp), float(max(depths)))
        #invert ticks values from positive to negative
        ax.annotate(
         AphaDict[ix],     # Text for the annotation
         xy = xy_coords,                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=FsizeA,                 # Font size for the annotation
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
    #Make the grid of the axis
     
   ##########################################################
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #fontsize of the yticks
    ax.tick_params(axis='both', labelsize=Figsize)
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
    ax.set_xlim(tstart, tend)
    #disable ticks on  x-axis
    ax.set_xticklabels([])
    #Free memory by  manually trigger garbage collection
    gc.collect()


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



def plot_spectrogram(ax, fig_object, time, tr, DatesList, Pressure= False, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
    #get sampling rate
    sampling_rate = tr.stats.sampling_rate
    #get the data
    data          = tr.data
    #The sampling rate is 1 second
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, DatesList)
#    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime' 
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
#    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    time          = mdates.date2num(time)
#    start_num     = min(time)
#    end_num       = max(time)
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
#   xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#   ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #Write the text on the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize= FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )

    #which axes object it should be near Calculate the colorbar position and size
    #######################
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    #cbar_x      = bbox.x1 + 0.09
    cbar_y      = bbox.y0
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #cbar_ax     = fig_object.add_axes([cbar_x , cbar_y -0.03, cbar_width, cbar_height])
    #set the colorbar annotation
    if(Pressure):
        #fig_object.colorbar(img, cax = cbar_ax, label = r'Power $(\mathrm{ dB \ ref\ Pa^2/Hz})$', format="%+2.f")
        fig_object.colorbar(img, cax = cbar_ax, label = r'$\mathrm{dB \ re\ 1Pa^2/Hz }$', format="%+2.f")
    else:
        #fig_object.colorbar(img, cax = cbar_ax, label = r'Power $(\mathrm{dB \ ref \ m^2s^{-4}/Hz})$', format="%+2.f")
        fig_object.colorbar(img, cax = cbar_ax, label = r'$\mathrm{dB \ re \ 1 m^2/s^{4}/Hz}$', format="%+2.f")
    #image_data = img.get_array()
    #Set font size for the label and ticks
    cbar_ax.tick_params(labelsize=Figsize)
    #Now set font size and padding through the colorbar's axis
    cbar_ax.yaxis.label.set_size(Figsize)        # Label font size
    cbar_ax.yaxis.labelpad = pad                  # Padding
#    PSD_MEANS  = np.mean(image_data, axis = 1)
    #################################
    #ENER_     = np.sum(np.square(spec))
#    ENER_      = np.mean(spec, axis =0)
    #ENER_      = np.square(spec)
#    ENRER_db   = 10 * np.log10(ENER_)
#    print(len(ENRER_db), len(PSD_MEANS))
#    file_name  = 'SEMIC_PDS_%s.dat'%(Date_of_Day)
    #Save the array into na file
#    np.savetxt(file_name , ENRER_db , delimiter="\n", fmt="%.4f")
    #Set label
#    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #set the y-axis as logarithm scale
#    ax.set_yscale('log')
#    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #disable ticks on  x-axis
    ax.set_xticklabels([])


def Func_PPSD(ax, fig_object, ListmseedFiles, fmxl, Pressure=False, **PARAMSDICT):
    #get the parameters  value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
    #Create an empty stream
    stream       = Stream()
    #Loop over the List of mseedFiles
    for fi,  mfile in zip(range(len(ListmseedFiles)), sorted(ListmseedFiles)):
        #add trace to a stream
        stream += read(mfile)
        #grap the Date
    #Merge the traces into one if they are continuous
    stream.merge(method=1)  # 'method=1' means interpolation if gaps are small
    #Load station metadata (StationXML file)
    #read the inventory file
    inventory   = read_inventory(fmxl)
    #Select first trace
    ######################################
    #Create PPSD object
    ppsd = PPSD(stream[0].stats, inventory)
    #if(Pressure):
    #    # Option 1: Without response correction (e.g., just use raw values)
    #    #ppsd = PPSD(stream[0].stats, metadata=None, special_handling="pressure")
    #    ppsd = PPSD(stream[0].stats, metadata=None)
    #else:
    #    ppsd = PPSD(stream[0].stats, inventory)
    ##Add entire 24-hour waveform data
    ppsd.add(stream)
    #Now we can delete the stream in order to free the memory
    del stream
    #clean the memory
    gc.collect()  # Force garbage collection
    #Original data
    psd_matrix  = np.array(ppsd.psd_values)       # shape: (n_times, n_periods)
    #Extract the periods from the ppsd
    periods     = np.array(ppsd.period_bin_centers)
    #Extract the time from the ppsd
    times       = [t.datetime for t in ppsd.times_processed]
    #Convert the times into the matplotlib
    time_nums   = mdates.date2num(times)
    #######################################################################
    #Transpose: shape (n_periods, n_times)
    Z_psd       = psd_matrix.T
    # --- Calculate extended time range: from 00:00 of day before to 00:00 of day after ---
    #first_time = times[0]
    #creates a new time at midnight (00:00) of the next day (24:00 boundary)
    #to Ensures your plot covers a full 24-hour period
    last_time   = times[-1]
    #00:00 of day before
    #start_time    = first_time.replace(hour=23, minute=59, second=0, microsecond=0) - timedelta(days=1)
    #start_time_num = mdates.date2num(start_time)
    #extended the last value of the times (last_time) to 24:00:00:00
    next_time   = last_time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(days=1)
    #end_time      = last_time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(days=1)
    #Convert the value to matplot time
    next_time_num = mdates.date2num(next_time) 
    #end_time_num  = mdates.date2num(end_time)
    #Adds a new axis (np.newaxis) to make it 2D so it can be horizontally stacked later
    #first_column  = Z_psd[:, 0][:, np.newaxis]    # shape: (n_periods, 1)
    last_column = Z_psd[:, -1][:, np.newaxis]             # shape: (n_periods, 1)
    #Horizontally appends the last column, repeating the last PSD value flat.
    #This is like saying “assume the last value continues until 24:00.”
    Z_extended  = np.hstack([Z_psd, last_column])          # shape: (n_periods, n_times + 1)
    #Z_extended = np.hstack([first_column, Z_psd, last_column])  # shape: (n_periods, n_times + 2)
    #Extends the time array to include the extra right edge at midnight.
    #This matches the new column added to Z.
    time_nums_extended  = np.append(time_nums, next_time_num)
    #Meshgrid
    #Creates 2D coordinate arrays for time and period.
    #T.shape and P.shape both match Z_extended.shape
    T, P                = np.meshgrid(time_nums_extended, periods)
    ###
    #T, P               = np.meshgrid(time_nums_extended, periods)

    #Draws the PSD spectrogram using the extended time and data.
    img  = ax.pcolormesh(T, P,  Z_extended,  shading='auto',
                              cmap=color)
    #shading='auto' avoids alignment issues and uses correct bin edges.

    #Set axis limits explicitly (to avoid showing extra data outside the range)
    #ax.set_xlim(time_nums_extended[0], time_nums_extended[-1])
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #Shift the first timestamp to midnight
    start_midnight = datetime.combine(times[0].date(), datetime.min.time())
    start_num      = mdates.date2num(start_midnight)
    #Keep your original end time
    end_num        = time_nums_extended[-1]
    ax.set_xlim(start_num, end_num)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)

    #set the y-axis as logarithm scale
    ax.set_yscale('log')
    #ax.set_ylabel("Period (s)")
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #Format the x-axis with dates
    ax.xaxis_date()
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height]) 
    #check if pressure
    if(Pressure):
        fig_object.colorbar(img, cax = cbar_ax, label = r'PPSD $(\mathrm{dB \ re\ 1Pa^2/Hz})$', format="%+2.f")
    else:
        fig_object.colorbar(img, cax = cbar_ax, label = r'PPSD $(\mathrm{dB \ re \ 1 m^2/s^{4}/Hz})$', format="%+2.f")
        ##################
    #Annotate the figure
    ax.annotate(
        AphaDict[ix],                   #Text for the annotation
        xy=  xy_coords,                 # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',        # Use axes fraction for positioning
        fontsize= FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
    #####################
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #########################
    #Set font size for the label and ticks
    cbar_ax.tick_params(labelsize=Figsize)
    #Now set font size and padding through the colorbar's axis
    cbar_ax.yaxis.label.set_size(Figsize)        # Label font size
    cbar_ax.yaxis.labelpad = pad                  # Padding
    #Label size
    #ax.tick_params(axis='both', labelsize=Figsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    #ax.set_ylim(1e-2, 1e3)
    #ax.set_ylabel([])
    #disable ticks on  x-axis
    ax.set_xticklabels([])




def plot_seismic_power(ax, time, tr, DatesList, Pressure=False, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
    #get sampling rate
    sampling_rate = tr.stats.sampling_rate
    #get the data
    #The sampling rate is 1 second
    #change the date_rng to numpy
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, DatesList)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime' 
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    ##############################33
    extent        = [start_num , end_num]
    # Calculate spectrogram
    #freqs, times, spec = spectrogram(data, fs=sampling_rate, nperseg=256, noverlap=128, scaling='density')
    freqs, times, spec = scipy.signal.spectrogram(tr.data, fs=1.0, nperseg=256, noverlap=128)
    ##########
    #get the coordinate of the point, where to write the text
    #which axes object it should be near Calculate the colorbar position and size
    #set the colorbar annotation
    #################################
    ENER_           = np.mean(spec, axis =0)
    #Call the interpolation
    time_new, ENER_ = inter1D(times, ENER_, time)
    #############################
    ENRER_db        = 10 * np.log10(ENER_)
    ##Plot Seismic power energy ###########
    ax.plot(time, ENRER_db, lw=8.0, linestyle ="-", color = 'grey', alpha =0.8)
    ax.plot(time, ENRER_db, lw=0.8, linestyle ="-", color = color, alpha =1.0)
    #Set label
    if(Pressure==True):
        #ax.set_ylabel(r'Power $(\mathrm{ dB \ ref\ Pa^2/Hz})$', labelpad = pad, fontsize = Figsize)
        ax.set_ylabel(r'$\mathrm{dB \ re\ 1Pa^2/Hz }$', labelpad = pad, fontsize = Figsize)
        #ax.set_ylabel('dB', labelpad = pad, fontsize = Figsize)
    else:
        ax.set_ylabel(r'$\mathrm{dB \ re \ 1 m^2/s^{4}/Hz}$', labelpad = pad, fontsize = Figsize)
    #ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #Set ylim of y-axis
    ax.set_ylim(min(ENRER_db), max(ENRER_db))
    #Set xlim of x-axis
    ax.set_xlim(min(time), max(time))
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         xy= xy_coords,                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=FsizeA,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #number of vertical lines for grid
    locations   = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Label size
    ax.tick_params(axis='both', labelsize=Figsize)
    #Make the grid of the axis
    ax.set_xlim(tstart, tend)
    #Set xlim of x-axis
    #ax.grid(visible = True, axis = "both", alpha = 0.7)
    ax.grid(visible = True, axis = "y", alpha = 0.5)
    ax.grid(visible = True, axis = "x", alpha = 0.4)
    #disable ticks on  x-axis
    ax.set_xticklabels([])

def Plot_Temp2D(ax, fig_object, time, DatesList, data2D, depths, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    vmind         = np.min(data2D, axis = 1)
    vmaxd         = np.max(data2D, axis = 1)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, DatesList)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime'
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    #Set some generic y-limits.
    y_lim_min     = np.min(depths)
    y_lim_max     = np.max(depths)
    #Set the extent for 2D plot
    extent        = [start_num , end_num,  y_lim_min, y_lim_max]
    #print(extent)
    #Make a 2D plot
    #im            = ax.imshow(data2D, cmap = color, vmin = 6.0 , vmax = 9.0, origin = 'lower', aspect = 'auto',
    im            = ax.imshow(data2D, cmap = color, origin = 'upper', aspect = 'auto',
                   interpolation ='spline16', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         xy = xy_coords,                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=FsizeA,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    #cbar_x      = bbox.x1 + 0.09
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #cbar_ax     = fig_object.add_axes([cbar_x , cbar_y -0.03, cbar_width, cbar_height])
    ########################
    #Set font size for the label and ticks
    cbar_ax.tick_params(labelsize=Figsize)
    #Add the colorbar on the figure
    #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,  format="%2.f")
    fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
    #wamba
    #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[6, 6.5,7 ,7.5,8.5,9])
    #fix the label position
    cbar_ax.yaxis.labelpad = pad
    #################################
    #Set label
    ax.set_ylabel('H (m)', labelpad = pad, fontsize = Figsize)
    #into a nice datetime string.
    #ax.xaxis_date()
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #if(ADCP_OPTION):
    #    #invert the y-axis
    #    ax.invert_yaxis()
    #    #invert ticks values from positive to negative
    #    ticks   = [-int(it) for it in ticks]
    #ax.invert_yaxis()
    #ticks   = [1 * int(it) for it in ticks]
    #ticks  =  np.arange(min(depths) , max(depths), 5)
    #ax.set_yticklabels(ticks)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    ax.set_ylim(float(min(depths)), float(max(depths)))
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #Make the grid
#    ax.yaxis.set_major_locator(mticker.MultipleLocator(base=20))
    ax.tick_params(axis='both', labelsize=Figsize)
    #disable ticks on  x-axis
    ax.set_xticklabels([])


def Extract_ADCPData(df,  Tlist,  param, average = None):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
    #Create an empty list to append the extracted data Check if the parameter==velocity
    if(check_if_string(Tlist)):
        Tlist =[ Tlist ]
    #Check if the desired velocity is the Horizontal velocity
    if(param=="velocity_up" or param=="velocity_down"):
        #Loop over the list of the time
        try:
            #Loop over the list of the time
            P_list = [np.nanmean(df.velds.U_mag.sel(time=ti), axis=0) for ti in Tlist]
            T_list = [df.velds.U_mag.sel(time=ti)['time'].values for ti in Tlist]
            ###############
            # Flatten lists and create DataFrame
            df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
            # Averaging duplicate timestamps
            df_extracted = df_extracted.groupby('time', as_index=False).mean()
            ####################################
            P            = df_extracted['value']
            T            = df_extracted['time'] 
            #################################
            del P_list ,T_list ,df_extracted 
            #clean the memory
            gc.collect()  # Force garbage collection
        except:
            print("The Data seems not to exist for the certain data in the list Time :  ")
            print(' , '.join(Tlist))
            exit()

    #Check if the desired velocity is the vertical velocity
    elif(param=="vertical_vel"):
        #Loop over the list of the time
        #Vertical component Z--vertical componet
#        P_list = [np.nanmean(df.velds.w.sel(time=ti), axis=0) for ti in Tlist]
#        T_list = [df.velds.w.sel(time=ti)['time'].values for ti in Tlist]
#        #North-Sud component N-S Component
        P_list = [np.nanmean(df.velds.v.sel(time=ti), axis=0) for ti in Tlist]
        T_list = [df.velds.v.sel(time=ti)['time'].values for ti in Tlist]
        #West-Est component  W-E Component
        #P_list = [np.nanmean(df.velds.u.sel(time=ti), axis=0) for ti in Tlist]
        #T_list = [df.velds.u.sel(time=ti)['time'].values for ti in Tlist]
        #Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list, T_list, df_extracted 
        #clean the memory
        gc.collect()  # Force garbage collection
    elif(param=="turbidity" or param=="avg_BS"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg  = np.mean(df.amp, axis = 0) 
        #Loop over the list of the time
        ########################
        P_list       = [np.nanmean(df_beam_avg.sel(time=ti), axis=0) for ti in Tlist]
        T_list       = [df_beam_avg.sel(time=ti)['time'].values for ti in Tlist]
        #Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list , T_list ,df_extracted 
        #clean the memory
        gc.collect()  # Force garbage collection
        ###########################
    elif(param=="dlnEcho"):
        #get the range form the data
        r           = df.range.data
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
#        Matrix2D    = np.asarray([df_beam_avg.sel(time=ti, range = r)  for ti in Tlist])
        Matrix2D   = xr.concat([df_beam_avg.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ####################################################
        T          = xr.concat([df_beam_avg.sel(time=ti)['time'] for ti in Tlist], dim="time")
#        #####################@
        DList       = np.copy(Matrix2D.values)
        #Free the memory for the large matrix
        del Matrix2D, df_beam_avg
        #clean the memory
        gc.collect()  # Force garbage collection
        #DList      = np.zeros(Matrix2D.values.shape, dtype = Matrix2D.values.dtype)
        for i in range(len(DList)):
            mean_val = np.mean(DList[i])
            DList[i] = [(item -mean_val)/mean_val for item in DList[i]]
        #DList = np.asarray(DList)
        if(average):
            ##%%%%%%%%%%%%% ##############
            #Return the Time and 1D perturbation in % (i.e. 100)
            return (T.values, np.mean(np.asarray(DList), axis=0) * 100) 
        ##############################
        #Return the Time and 2D perturbation in % (i.e. 100)
        return (T.values, np.asarray(DList) * 100)
        #################################################### 
    elif(param=="velocity2D_UP" or param=="velocity2D_DOWN"):
        #get the range form the data
        r         = df.range.data
        #get the Horizontal velocity
        U         = df.velds.U_mag
        #extract the 2D velocity
        Matrix2D  = xr.concat([U.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ###################################################
        T         = xr.concat([U.sel(time=ti)['time'] for ti in Tlist], dim="time")
        #####################@
        #Free the memory for the large matrix
        del U, r
        #clean the memory
        gc.collect()  # Force garbage collection
        #return the values
        return (T.values, Matrix2D.values)
    ##%%%%%%%% Extract the current direction %%%%%%%%%%%%%%%%%%
    elif(param=="veldir2D_UP" or  param=="veldir2D_DOWN"): 
        #get the range form the data
        r        = df.range.data
        #Get the velocity Direction in 2D
        U        = df.velds.U_dir
        #extract the 2D velocity
        Matrix2D = xr.concat([U.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ####################################################
        T        = xr.concat([U.sel(time=ti)['time'] for ti in Tlist], dim="time")
        #####################@
        #Free the memory for the large matrix
        del U, r
        #clean the memory
        gc.collect()  # Force garbage collection
        #return the values
        return (T.values, Matrix2D.values)
    ################################################
    else:
        #Loop over the  list of the time an extract the desire parameter
        P_list       = [df[param].sel(time=ti) for ti in Tlist]
        T_list       = [df[param].sel(time=ti)['time'].values for ti in Tlist]
        # Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()  # Averaging duplicate timestamps
        #################
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list, T_list, df_extracted 
    #Return the  value corresponding to the date of your choice
    return (T, P)

#####%%%%%%%% Extract the Temperature %%%%%%%%%%%%%%%%%%%%%%%%%%%
def Extract_RBR_SOLO_TEMP(DataDict,  Tlist, AVERAGE=False):
    #The data should be in the dictionary, where the keys are the depths
    DATA2D     = []
    DEPTHS     = []
    #get the time you need to plot
    TIMES      = set([time for depth in DataDict.keys() for time in DataDict[depth] if(time.split()[0] in Tlist)])
    #####################################
    TIMES      = sorted(np.asarray(list(TIMES)))
    #loop over the depths and the time to form a matrix
    for  depth in DataDict.keys():
        #get the temperature at a specific depth
        temp   = [DataDict[depth][time] for time in TIMES]
        #print(temp)
        DATA2D.append(temp)
        DEPTHS.append(float(depth))
    #Trans formed the List into# 2D array
    DATA2D     = np.asarray(DATA2D, dtype=float)
    #Convert the TIMES into pan#das time
    TIMES      = pd.to_datetime(TIMES)
    #Check if you need everage #at the depth
    if(AVERAGE):
        #data   = np.mean(DATA2D, axis = 0)
        data   = np.nanmean(DATA2D, axis = 0)
        return (TIMES, data)
    else:
        return(TIMES,DEPTHS, DATA2D)



##############plot histogram ######################
def Plot_hist(ax, time, data, param, DatesList,  **PARAMSDICT):
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    nbins             = hist_nbin(data)
    #print(nbins)
    #ax.bar(time, data, width = 0.01, align ='edge')
    #compute the width as function of the data
    #xwd = 1/data.size + (1./data.size) * 0.89
    if(len(set(data)) > 2):
        xwd = 1/data.size + (1./data.size) * 0.5
        xwd = 0.0067
        #print(xwd)
        #exit()
        #xwd = 1/data.size + (1./data.size) * 0.1
        ax.bar(time, data, width = xwd , align ='edge', color = color, lw= 100.0, label = param.title())
        #ax.bar(time, data, width = xwd , align ='edge', color = color, lw= 50.0, label = param.title())
    else:
        ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label = param.title())
    #ax.legend(loc="upper left", prop = {'size': 6})
    #Annotate the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=   xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize= FsizeA,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    ###################
    #get the start and the endtime
    tstart, tend  =  start_endtime(time, DatesList)
    ###################################
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #ax.yticks(fontsize = fsize)
    ax.tick_params(axis='y', labelsize= Figsize)
    ax.set_xlim(tstart, tend)
    #Set label
    ax.set_ylabel(ylabel,  fontsize = Figsize)
    ################################
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
with open("confsevdays_stck.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Grab the Dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFiles          =   Fig_params['mseedFiles']
###########################################
#STARTDATE      =   Fig_params['STARTDATE']
#Seimic mseed file
Response_File       = Fig_params['Response_File']
#%%%%%%%%%%Pressure%%%%%%%%%%%%%%%%%%%%%%%%
#Grab the path of the data
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
##########################################
dlnEcho_DOWN        = Fig_params['dlnEcho_DOWN']
dlnEcho_UP          = Fig_params['dlnEcho_UP']
#ADCP parameter to plot
param_adpc          = Fig_params['param_adpc']
velocity_up         = Fig_params['velocity_up']
velocity_down       = Fig_params['velocity_down']
#2D velocity
velocity2D_UP       = Fig_params['velocity2D_UP']
velocity2D_DOWN     = Fig_params['velocity2D_DOWN']
#Parameter to plot Vertical  Velocity, this actually the Vertical velocty given the componant w
vertical_vel        = Fig_params['vertical_vel']
#Do want to plot water Pressure?
pressure_up         = Fig_params['pressure_up'] 
pressure_down       = Fig_params['pressure_down'] 
#Plot 2D current direction
veldir2D_UP         = Fig_params['veldir2D_UP']
veldir2D_DOWN       = Fig_params['veldir2D_DOWN'] 
#Discharge file name
FILE_DISCH         = Fig_params['FILE_DISCH']
##Lake Level File
FILE_LAKE          = Fig_params['FILE_LAKE']
#Grab the Meteo data file
FILE_METEO     = Fig_params['FILE_METEO']
#Grab the temperature from BRB-SOLO-3
FILE_TEMPEA     = Fig_params['FILE_TEMPEA']
FILE_TEMPEB     = Fig_params['FILE_TEMPEB']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']
#############################################
#Check the user want to make a plot bigger?
#plot_bigger    = Fig_params['plot_bigger']
wind_speed        = Fig_params['wind_speed']
#############################################
#Plot the precipitation
precipitation  = Fig_params['precipitation']
#Atmospheric pressure
Atm_pressure_QFE  = Fig_params['Atm_pressure_QFE']
Atm_pressure_QFF  = Fig_params['Atm_pressure_QFF']
#get the Air temperature parameter
TemperatureAir    = Fig_params['TemperatureAir']
Atm_pressure_Temp = Fig_params['Atm_pressure_Temp']
#########################################
Lake_Level        = Fig_params['Lake_Level']
#load the discharge
discharge         = Fig_params['discharge']
# Plot the seismic Waveform 
waveform          = Fig_params['waveform']
#get the envelope
envelope          = Fig_params['envelope']
#######################
seismic_power     = Fig_params['seismic_power']
#######################
wind_seismic_power = Fig_params['wind_seismic_power']
#########################################
seismicpower_dlnEcho  = Fig_params['seismicpower_dlnEcho']
#########################
pressure_waveform = Fig_params['pressure_waveform']
######## plot the spectrogram ###
spectrogram       = Fig_params['spectrogram']
#Plot seismic spectrogram probabilitic spectral density
PPSD_spec         = Fig_params['PPSD_spec']
#ADCP temperature
temperature_up    = Fig_params['temperature_up']
temperature_down  =  Fig_params['temperature_down']
#get the Thermistors temperature
Temperature2D_MA  = Fig_params['Temperature2D_MA']
Temperature2D_MB  = Fig_params['Temperature2D_MB']
Temperature_AVG   = Fig_params['Temperature_AVG']
#Grab the bandpass
bandpass          = Fig_params['bandpass']
#bandpass full spectrum
bandpass_spec     = Fig_params['bandpass_spec']
#Plot in a window, set this true if you're plottion a window in a one day
Windowing         = Fig_params['Windowing']
#Set the number of subfigures to Zero
##change the values on the number of figure to plot by re-signing the len of PARAMSDICT
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
#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size,constrained_layout=True)
#fig, axs  = plt.subplots(nfigs, 1, sharex = True, figsize = fig_size,constrained_layout=True,gridspec_kw={'hspace': 0.02})
#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
########################################################################################@
#Read the Turbidity  file
#d_disch= pd.read_fwf(FILE_DISCH, delimiter='\t')
##Add the of DateTime
#d_disch= Add_Column(d_disch)
#Load the Discahge file
d_disch= pd.read_csv(FILE_DISCH, skiprows=0, sep =",",encoding="latin1")
#print(d_disch) 
#read the Meteo database
d_mteo = pd.read_csv(FILE_METEO, delimiter ='\t')
#Add the of DateTime
d_mteo = Add_Column(d_mteo)

#load the Lake Level Data
dfLake = pd.read_csv(FILE_LAKE, skiprows=0, sep =",",encoding="latin1")
#Read the discharge data file
#Load the Discahge file
d_disch= pd.read_csv(FILE_DISCH, skiprows=0, sep =",",encoding="latin1")
#Grab 2D Temperature
#get the temperature data from RBR-SOLO-3
# Opening JSON file
ftmpA  = open(FILE_TEMPEA)
ftmpB  = open(FILE_TEMPEB)
# returns JSON object as
# a dictionary
data_tempA = json.load(ftmpA)
data_tempB = json.load(ftmpB)

####### Read the ADCP FILES##############
try:
    df_down = dlfn.read(ADCP_FILE_NAME_DOWN)
    df_up   = dlfn.read(ADCP_FILE_NAME_UP)
except:
    print("print provide the ADCP files:  %s and %s  "%(ADCP_FILE_NAME_UP, ADCP_FILE_NAME_DOWN))
#plot_discharge(ax1,hours_of_day, data_day, starttime, Remove_frame)

#open an empty Date set
DictDate     = set()
#Create an empty stream
stream       = Stream()
#stream of spectrogram
stream_spec = Stream()
#set the pressure variable here to be same nature a pressure_waveform
Pressure     = pressure_waveform
#loop over the mseedFiles List
for ix,  mseedFile in zip(range(len(mseedFiles)), sorted(mseedFiles)):
    #plot_bigger = False
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
    #set the parameter
    #get the seismogram from mseedFile
    time, tr, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass, Pressure)
    #add trace to a stream
    stream += tr
    #check you want to plot the spectrum
    if(spectrogram):
        time, tr_spec, Title = read_mseed(DAYDATE, mseedFile, Response_File, bandpass_spec, Pressure)
        stream_spec += tr_spec
    #grap the Date
    DictDate.add(DAYDATE)
#print(stream)
#Merge the traces into one if they are continuous
if(len(mseedFiles) > 1):
    stream.merge(method=1)  # 'method=1' means interpolation if gaps are small
    #try this
######### Merge the spectrum #################
if(spectrogram and len(mseedFiles) > 1):
    stream_spec.merge(method=1)
#Optional: Check for any gaps or overlaps
#Transform the Date set into the list
DictDate = list(sorted(DictDate))
#print(stream.get_gaps())
#find the time for the correspond trace
date_rng = pd.date_range(start= DictDate[0], freq='1s', periods= stream[0].data.size)
#change the date_rng to numpy
time_tr  = date_rng.to_numpy()
#compute the seismic power
Dict_idx= set()
#plot Waveform
if(waveform):
    #get the Dischage data for the corresponding time
    param         = "waveform"
    #get the Dischage data for the corresponding time
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    if(ix not in Dict_idx):
        tle_obj = axs[ix].set_title(station, loc='center', fontsize = Figsize, pad=None)
        #tle_obj = axs[ix].set_title(station, loc='center',fontsize =8, pad=-50)
        #tle_obj.set_zorder(10)
        #tle_obj.set_zorder(5)
        Dict_idx.add(ix)
    #Plot the figure by calling the plotting function, plot_twinx
    #Plot_fig(axs[ix], time_tr, stream[0].data * 1e+6, DictDate, param, Plot_small= True, **PARAMSDICT)
    tr_filt = stream[0].copy()
    #tr_filt.filter('lowpass', freq=0.05, corners=2, zerophase=True)
    tr_filt.filter('lowpass', freq=0.03, corners=2, zerophase=True)
    #tr_filt.filter('lowpass', freq=0.03, corners=4, zerophase=False)
    #tr_filt.filter('highpass', freq=0.1, corners=2, zerophase=True)
    Plot_fig(axs[ix], time_tr, tr_filt.data * 1e+6, DictDate, param, Plot_small= True, **PARAMSDICT)
    ###############################################
    #Plot_fig(axs[ix], time_tr, stream[0].data * 1e+6 / max(stream[0].data * 1e+6), DictDate, param, Plot_small= True, **PARAMSDICT)
    #set the limit
    #axs[ix].set_ylim(-0.5, 0.5)
#    axs[ix].set_ylim(-15, 15)
    #bring axs in front
    #Bring x-axis to the front
#    axs[ix].set_axisbelow(False)  # Prevents the grid and axes from staying behind plots
#    axs[ix].spines['bottom'].set_zorder(15)  # Ensures x-axis is drawn on top
    # Bring second subplot (ax2) to front using zorder
    #for spine in axs[ix].spines.values():
    #    spine.set_zorder(100)
    ###########################
    #axs[ix].set_zorder(50)
    #axs[ix].patch.set_visible(False)  # To avoid hiding the first plot's background


if(pressure_waveform):
    #get the Dischage data for the corresponding time
    param         = "pressure_waveform"
    #get the Dischage data for the corresponding time
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    if(ix not in Dict_idx):
        tle_obj = axs[ix].set_title(station, loc='center', fontsize = Figsize, pad=None)
        tle_obj.set_zorder(10)
        Dict_idx.add(ix)
    #Plot the figure by calling the plotting function, plot_twinx
    #Plot_fig(axs[ix], time, trace.data, list(DictDate), param, **PARAMSDICT)
    Plot_fig(axs[ix], time_tr, stream[0].data, DictDate, param, Plot_small= True, **PARAMSDICT)
    #Zomm on the plot
#    axs[ix].set_ylim(-15, 15)
#    axs[ix].set_axisbelow(False)  # Prevents the grid and axes from staying behind plots
#    axs[ix].spines['bottom'].set_zorder(15)  # Ensures x-axis is drawn on top


if(envelope):
    #set the corresponding parameter
    param         = "envelope"
    #get the Dischage data for the corresponding time
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    ##Plot the envelope
    PlotEnvelop(axs[ix], time_tr, stream[0], DictDate, param,  Pressure, linewidth=1.2, **PARAMSDICT)

#plot the spectrogram
if(spectrogram):
    #set the parameter
    param            = "spectrogram"
    #get the seismogram from mseedFile
    #if(pressure_waveform):
    #    Pressure=True
    #else:
    #    Pressure=False
    ##get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #plot_spectrogram(axs[ix], fig, trace.data, list(DictDate), Pressure, **PARAMSDICT)
    plot_spectrogram(axs[ix], fig, time_tr, stream_spec[0], DictDate, Pressure, **PARAMSDICT)
    #zoom in the spectrogram
#    axs[ix].set_ylim(0, 0.25)

if(PPSD_spec):
    #set the parameter
    param            = "PPSD_spec"
    ##get the seismogram from mseedFile
    ##get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Func_PPSD(axs[ix], fig, mseedFiles, Response_File, Pressure, **PARAMSDICT) 


if(seismic_power):
    #set the parameter
    param            = "seismic_power"
    #get the seismogram from mseedFile
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the seismic energy
    plot_seismic_power(axs[ix], time_tr, stream[0], DictDate, Pressure, **PARAMSDICT)
#DictDate = ["2023-10-22", "2023-10-23"]
#Plot the wind and the seimsic power together
if(wind_seismic_power):
    #set the parameter
    param              = "wind_seismic_power"
    #get the wind speed data for the corresponding time
    time, data         = Extract_database(d_mteo,  DictDate, "wind_speed")
    #get the seismic energy calculation
    freqs, times, spec = scipy.signal.spectrogram(stream[0].data, fs=1.0, nperseg=256, noverlap=128)
    #compute the energy
    ENER_              = np.mean(spec, axis =0)
    #Call the interpolation
    time_new, ENER_    = inter1D(times, ENER_, time)
    #Compute the Energy in decibel
    ENRER_db           = 10 * np.log10(ENER_)
    ##Plot Seismic power energy ###########
    #Sort and drop duplicates if needed
    df = pd.DataFrame({'time': time})
    df = df.groupby('time', as_index=False).mean()
    #########################
    #get the parameter, the index of the corresponding axis, and the color
    label_ , color, ix   = PARAMSDICT[param]
    #Plot the wind speed
    Plot_fig(axs[ix], df['time'], data, DictDate, param, **PARAMSDICT)
    #set the color of the y-label same as the curve color
    axs[ix].tick_params(axis='y', colors=color)
    axs[ix].set_ylabel(axs[ix].get_ylabel(), color=color)
    #Create a second y-axis sharing the same x-axis
    ax2     = axs[ix].twinx()
    ##############
    #set a new color
    power_color = 'k'
    #Copy the dictionary
    new_dict = PARAMSDICT.copy()
    if(pressure_waveform):
        ylabel_new = r'$\mathrm{dB \ re\ 1Pa^2/Hz }$'
    else:
        ylabel_new = r'$\mathrm{dB \ re \ 1 m^2/s^{4}/Hz}$'
    ###########################
    #Make changes in the dictionary
    new_dict[param][0] = ylabel_new 
    new_dict[param][1] = power_color 
    #Sort and drop duplicates if needed
    ### %%%%%%%%%%%%%%%%%% ##########
    #Do not change this, only if you understand what you are doing
    Plot_fig(ax2, df['time'], ENRER_db, DictDate, param, **new_dict)
    #### %%%%%%%%%%%%%%%%%%%%
    #set the label
    ax2.set_ylabel(ylabel_new, color=power_color, labelpad = pad, fontsize=Figsize)
    #remove the grid on the axis
    ax2.grid(axis='both', visible=False)
    #set the color of the  ticks and label
    ax2.tick_params(axis='y', colors=power_color)
    #format the y-axis to be intergers
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:d}'.format(int(val))))
    #get the start and the endtime
    #disable axis
    ax2.set_xticklabels([])

if(Atm_pressure_QFE):
    #get the Dischage data for the corresponding time
    param         = "Atm_pressure_QFE"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  DictDate, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

if(Atm_pressure_QFF):
    #get the pressure atmospheric data for the corresponding time
    param         = "Atm_pressure_QFF"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  DictDate, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

if(TemperatureAir):
    #get the Air Temperature data for the corresponding time
    param         = "TemperatureAir"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  DictDate, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)


if(Atm_pressure_Temp):
    #get the atmospheric pressure and Air Temperature data for the corresponding time
    param           = "Atm_pressure_Temp"
    #get the atmospheric pressure data for the corresponding time
    timep, datap    = Extract_database(d_mteo,  DictDate, "Atm_pressure_QFF")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix   = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], timep, datap, DictDate, param, **PARAMSDICT)
    #####################################################################
    #Create a second y-axis sharing the same x-axis
    ax2             = axs[ix].twinx()
    ##############
    #set a new color
    bt_color        = 'r'
    #ylabel_echo = 'd'+r'$\ln$'+'(Echo)'
    ylabel_Temp     = "Air Temp (°C)"
    ###########################
    #Copy the dictionary
    new_dict = PARAMSDICT.copy()
    #Make changes in the dictionary
    new_dict[param][0] = ylabel_Temp 
    new_dict[param][1] = bt_color 
    ### %%%%%%%%%%%%%%%%%% ##########
    #get the Air Temperature data for the corresponding time
    timeT, dataT       = Extract_database(d_mteo,  DictDate, "TemperatureAir")
    #Do not change this, only if you understand what you are doing
    Plot_fig(ax2, timeT, dataT, DictDate, param, **new_dict)
    #### %%%%%%%%%%%%%%%%%%%%
    #set the label
    ax2.set_ylabel(ylabel_Temp, color=bt_color, labelpad = pad, fontsize=Figsize)
    #set ticks color
    ax2.tick_params(axis='y', colors=bt_color)
    #remove the grid on the axis
    ax2.grid(axis='both', visible=False)
    #get the start and the endtime
    #disable axis
    ax2.set_xticklabels([])


if(wind_speed):
    #get the Dischage data for the corresponding time
    param         = "wind_speed"
    #get the wind speed data for the corresponding time
    time, data    = Extract_database(d_mteo,  DictDate, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
    #Plot_fig(axs[ix], time, detrend(data), list(DictDate), param, **PARAMSDICT)

#Plot the precipitation
if(precipitation):
    #set the parameter
    param         = "precipitation"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  DictDate, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_hist(axs[ix], time, data, param, DictDate, **PARAMSDICT)


##########Check if the user want to plot the sea level
if(Lake_Level):
    param         = "Lake_Level"
    #get the Dischage data for the corresponding time
    time, data    = ExtractLakeLevel(dfLake,  DictDate)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    #Plot_fig(axs[ix], time, data, param, plot_bigger,  **PARAMSDICT)
    Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

##Plot the Discharge
if(discharge):
        #set the parameter
        param          = "discharge"
        #get the Dischage data for the corresponding time
        #time, data    = Extract_database(d_disch,  date_list, 'Discharge')
        time, data     = ExtractLakeLevel(d_disch,  DictDate)
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix  = PARAMSDICT[param]
        #Plot the figure by calling the plotting function, plot_twinx
        Plot_fig(axs[ix], time, data, DictDate,  param, **PARAMSDICT)

if(Temperature2D_MA):
   param            = "Temperature2D_MA"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, depths, data2D       = Extract_RBR_SOLO_TEMP(data_tempA,  DictDate, AVERAGE=False)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix            = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_Temp2D(axs[ix],fig, time, DictDate, data2D, depths,  **PARAMSDICT)
   #Plot_Temp2D(axs,fig, time, data2D, depths, **PARAMSDICT)

if(Temperature2D_MB):
   param            = "Temperature2D_MB"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, depths, data2D       = Extract_RBR_SOLO_TEMP(data_tempB,  DictDate, AVERAGE=False)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix            = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_Temp2D(axs[ix],fig, time, DictDate, data2D, depths,  **PARAMSDICT)

if(Temperature_AVG):
   param            = "Temperature_AVG"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time,  data       = Extract_RBR_SOLO_TEMP(data_tempB,  DictDate, AVERAGE=True)
   #time,  data       = Extract_RBR_SOLO_TEMP(data_tempA,  DictDate, AVERAGE=True)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix            = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, DictDate,  param, **PARAMSDICT)

##Doubleched if we need to plot the Temperature
if(temperature_up):
   #set the parameter
   param        = "temperature_up"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_ADCPData(df_up , DictDate, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix= PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, DictDate,  param, **PARAMSDICT)

if(temperature_down):
   #set the parameter
   param        = "temperature_down"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_ADCPData(df_down , DictDate, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix= PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

##Doubleched if we need to plot the Velocity
if(velocity_up):
   #set the parameter
   param      = "velocity_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_ADCPData(df_up , DictDate, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
   #Plot the figure by calling the plotting function
if(velocity_down):
   #set the parameter
   param      = "velocity_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_ADCPData(df_down , DictDate, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
   #Plot the figure by calling the plotting function
if(vertical_vel):
       #set the parameter
       param      = "vertical_vel"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_ADCPData(df_up , DictDate, param)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix = PARAMSDICT[param]
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
##Plot the 2D velocity
if(velocity2D_UP):
    #set the parameter
    param      = "velocity2D_UP"
    #get the Turbidity data for the corresponding time
    time, data2D = Extract_ADCPData(df_up, DictDate, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix], fig, df_up, time, data2D, DictDate, ADCP_UP_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="top" , **PARAMSDICT)

##Plot the 2D velocity
if(velocity2D_DOWN):
    #set the parameter
    param      = "velocity2D_DOWN"
    #get the Turbidity data for the corresponding time
    time, data2D = Extract_ADCPData(df_down, DictDate, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix], fig, df_down, time, data2D, DictDate, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=True, rframe="top" , **PARAMSDICT)

##Plot the 2D Current Direction
if(veldir2D_UP):
    #set the parameter
    param      = "veldir2D_UP"
    #get the Turbidity data for the corresponding time
    time, data2D = Extract_ADCPData(df_up, DictDate, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix], fig, df_up, time, data2D, DictDate, ADCP_UP_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="top" , **PARAMSDICT)

##Plot the 2D Current Direction
if(veldir2D_DOWN):
    #set the parameter
    param      = "veldir2D_DOWN"
    #get the Turbidity data for the corresponding time
    time, data2D = Extract_ADCPData(df_up, DictDate, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix], fig, df_down, time, data2D, DictDate, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=True, rframe="top" , **PARAMSDICT)

if(pressure_up):
       #set the parameter
       param      = "pressure"
       #grab the data
       time, data = Extract_ADCPData(df_up , DictDate, param)
       #modify the parameter
       param      = param+"_up"
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix = PARAMSDICT[param]
       #Plot the figure by calling the plotting function
       #convert decibar to bar; 1dbar = 0.1 bar
       #data  = 0.1 * data
       #convert decibar to pascal; 1dbar = 10000 Pa
#       data  = 10000 * data
#       #make a new dictionary
       #new_dict = PARAMSDICT.copy()
       PARAMSDICT[param][0]= 'dbar'
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

if(pressure_down):
       #set the parameter
       param      = "pressure"
       #Grab the data
       time, data = Extract_ADCPData(df_down , DictDate, param)
       #modify the parameter
       param      = param+"_down"
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix = PARAMSDICT[param]
       #Plot the figure by calling the plotting function
       #convert decibar to bar; 1dbar = 0.1 bar
       data  = 0.1 * data
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)

###################################
if(dlnEcho_UP):
    #set the parameter
    param      = "dlnEcho"
    #Grab the data
    time, Matrix2D = Extract_ADCPData(df_up , DictDate, param)
    #################
    #change the paramter
    param      = "dlnEcho_%s"%("UP")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #plot the 2D Backscatter
    #plot the upward looking ADCP
    #Plot_fig2D(axs[ix], fig, df_up, time, DictDate, Matrix2D, ADCP_UP_HEIGTH_FROM_LAKEBED, param, rframe= None, **PARAMSDICT)

    Plot_fig2D(axs[ix], fig, df_up, time, Matrix2D, DictDate, ADCP_UP_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="top" , **PARAMSDICT)

if(dlnEcho_DOWN):
    #set the parameter
    param      = "dlnEcho"
    #grab the data
    time, Matrix2D = Extract_ADCPData(df_down , DictDate, param)
    #changed the paramter
    param      = "dlnEcho_%s"%("DOWN")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #plot the 2D Backscatter
    #plot the upward looking ADCP
    #Plot_fig2D(axs[ix], fig, df_down, time, DictDate, Matrix2D,ADCP_DOWN_HEIGTH_FROM_LAKEBED, param, rframe= None, **PARAMSDICT)

    Plot_fig2D(axs[ix], fig, df_down, time, Matrix2D, DictDate, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=True, rframe="top" , **PARAMSDICT)
#set the x-label
if(seismicpower_dlnEcho):
    param      = "seismicpower_dlnEcho"
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #grab the data
    time_BS, dlnBS_down = Extract_ADCPData(df_down , DictDate, "dlnEcho", average=True)
    #Convert pandas datetime to matplotlib date format
    tadd       = mdates.date2num(time_BS)
    #Compute the spectrogram
    freqs, times, spec = scipy.signal.spectrogram(stream[0].data, fs=1.0, nperseg=256, noverlap=128)
    ##########
    #plot the backscatter perturbation
    #Copy the dictionary
    new_dict = PARAMSDICT.copy()
    if(pressure_waveform):
        PARAMSDICT[param][0]= 'Pa²/Hz'
    #set the colorbar annotation
    #################################
    ENER_           = np.mean(spec, axis =0)
    #Call the interpolation (downsampling) to fit the ADCP data
    #############################
    ENRER_db        = 10 * np.log10(ENER_)
    #interpolate the energy
    eTime, eData    =  inter1DD(times,  ENRER_db, tadd)
    ##Plot Seismic power energy ###########
    #Sort and drop duplicates if needed
    #df = pd.DataFrame({'time': time_BS, 'value': eData})
    df = pd.DataFrame({'time': time_BS})
    df = df.groupby('time', as_index=False).mean()
    ##########################
    #Do not change this, only if you understand what you are doing
    Plot_fig(axs[ix], df['time'], eData, DictDate, param,  **PARAMSDICT)
    ###############
    axs[ix].yaxis.set_major_formatter(FormatStrFormatter('%d'))
    #Create a second y-axis sharing the same x-axis
    ax2     = axs[ix].twinx()
    ##############
    #set a new color
    bs_color = 'r'
    ylabel_echo = 'd'+r'$\ln$'+'(Echo)'
    ###########################
    #Make changes in the dictionary
    new_dict[param][0] = ylabel_echo    
    new_dict[param][1] = bs_color 
    #Sort and drop duplicates if needed
    ### %%%%%%%%%%%%%%%%%% ##########
    #Do not change this, only if you understand what you are doing
    Plot_fig(ax2, df['time'], dlnBS_down, DictDate, param, **new_dict)
    #### %%%%%%%%%%%%%%%%%%%%
    #set the label
    ax2.set_ylabel(ylabel_echo, color=bs_color, labelpad = pad, fontsize=Figsize)
    #set ticks color
    ax2.tick_params(axis='y', colors=bs_color)
    #remove the grid on the axis
    ax2.grid(axis='both', visible=False)
    #get the start and the endtime
    #disable axis
    ax2.set_xticklabels([])

#Write text on the x-axis at specific_time 
#axs[nfigs-1].text(time_write_1, 2.0, 'Wind onset', rotation=0, color='k')
#axs[nfigs-1].text(time_write_2, 2.0, 'Seismic signal onset', rotation=0, color='k')


#window selection


#Formatted the x-axis
if(len(DictDate) ==1):
    #figname   = "Waveform_%s.png"%(list(DictDate)[0])
    if(pressure_waveform):
        figname   = "Pressure_Waveform_%s.pdf"%(list(DictDate)[0])
    else:
        figname   = "Waveform_%s.pdf"%(list(DictDate)[0])
    ###%%%%%%%%%% set the date on the x-axis %%%%%%%%%%%%%%%%%
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #### Adjust font size of the x-axis labels
    axs[-1].tick_params(axis="x", labelsize=Figsize)      
    #write on the last x-axis
    #label the x-label
    axs[-1].set_xlabel('Time (hour:minute) on %s'%(list(DictDate)[0]), labelpad =12, fontsize = Figsize)
else:
    #figname   = "Waveform_%s_%s.png"%(list(DictDate)[0],list(DictDate)[-1])
    figname   = "Waveform_%s_%s.pdf"%(list(DictDate)[0],list(DictDate)[-1])
    ###%%%%%%%%%% set the date on the x-axis %%%%%%%%%%%%%%%%%
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    axs[-1].tick_params(axis='x', rotation=45)



##############################################
plt.xticks(fontsize= 18)
########################################
# Bring the window to the front (works with Qt5Agg or TkAgg backends)
#plt.subplots_adjust(hspace = 0.08)
plt.subplots_adjust(hspace = fig_space)
#Align ylabel of the figure
fig.align_ylabels(axs)
#Save the figure
#fig.tight_layout()
#fig.canvas.manager.window.raise_()
#Save the figure
fig.savefig(figname, bbox_inches = 'tight')
