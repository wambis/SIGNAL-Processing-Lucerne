#!/usr/bin/env python

#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
import obspy
from obspy import read, read_inventory
from obspy.signal import PPSD
from glob import glob
from scipy import stats
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from math import log10, floor, ceil
import numpy as np
import pandas as pd
from tabulate import tabulate
#for interpolation
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import interpolate
import xlrd
#################################################
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
##############################################
import matplotlib.dates as dt
import xarray as xr
import pandas as pd
from mhkit import dolfyn
from mhkit import dolfyn as dlfn
from datetime import datetime, timedelta
#################
#from scipy.signal import spectrogram
import scipy
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
import matplotlib.dates as mdates
from mhkit.dolfyn.adp import api
from pandas.core.common import flatten
##################################
from obspy.core import UTCDateTime




#AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
#            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
#            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}

#AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
#            5: 'f', 6: 'g', 7:'h', 8:'i', 9: 'j',10:'k',
#            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}

AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            5: 'f', 6: 'f', 7:'h', 8:'g', 9: 'j',10:'k',
            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}
#size for anotation
f_size = 13
Figsize = 12
pad = 10
#############
xy_coords=(0.0089, 0.85)                   # Coordinates (relative to axes) for the annotation
############################################

def inter1D(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations


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

#def Scatt_Correct(df, Profile_range, angle = 20.0):
def BSC_Correct(df, Tlist, alpha =0.6, nrange=25, MEAN = True):
    #get the parameter
    beam_angle = df.velds.attrs['beam_angle']
    #get the Profile range
    PR         = df.range.data[: nrange]
    #Calculate the range along the beam
    R          = PR * np.cos(beam_angle)
    #########################################
    #amp_new    = df.amp * 0.43 + 20 * np.log10(R) + 2 *alpha * R
    #alpha is water absorption
    #alpha      = 0.6
    #alpha      = 2.4  #Good value for both ADCP (UP and down)
    ####################################
    #make the average on all the 4 beams, NB the key word here is amp
    df_beam_avg = np.mean(df.amp, axis = 0) 
    #############################
    if(MEAN):
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
        #Correct the Bascatter
        PC = P * 0.43 + 20 * np.log10(R) + 2 *alpha * R
    else:
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
        #Get the range, altitude
        #r           =  df_beam_avg.range
        r           =  df_beam_avg.range.data[:nrange]
        #Loop over the list of the time
        Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist], axis =1 )
        #time        = np.concatenate([df.amp.sel(time = ti) for ti in Tlist])
        time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
        M2D_NEW     = []
        #for ir, slide in zip(np.arange(len(Matrix2D)), Matrix2D):
        for r_i, slide in zip(R, Matrix2D):
            #print(i, slide)
            slide_new = slide *0.43 +  20 * np.log10(r_i) + 2 *alpha * r_i
            #print(slide_new)
            M2D_NEW.append(slide_new)
            #print(type(slide_new))
        #Correct the Bascatter
        #PC = Matrix2D * 0.43 + 20 * np.log10(R) + 2 *alpha * R
        PC  = np.asarray(M2D_NEW)
    return time, PC

def xy_point(xmin, xmax, ymin, ymax):
    xp         = (xmin + xmax)/2
    #xp         = xmin + 0.28 * (xmax - xmin)
    #yp         = ymin + 0.8 * (ymax -ymin)
    yp         = ymin + 0.85 * (ymax -ymin)
    return(xp, yp)

#Define the Plot for the Envelop

def PlotEnvelop(ax, Date_of_Day, param, tr_filt, linewidth=1.2,  **PARAMSDICT):
    #get the sampling rate
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    sampling_rate = tr_filt.stats.sampling_rate

    #get the data from the filtered trace
    #data   = tr_filt.data

    #tr_filt.filter('bandpass', freqmin=1.0/200, freqmax=1.0/167, corners=2, zerophase=True) #Good for 2023-10-23
    tr_filt.filter('bandpass', freqmin=1.0/200, freqmax= 1.0/158, corners=2, zerophase=True) #Good for 2023-10-23
    #tr_filt.filter('bandpass', freqmin=0.006, freqmax= 0.015, corners=2, zerophase=True) #Good for 2023-09-02
#    tr_filt.filter('bandpass', freqmin=0.004, freqmax= 0.01, corners=2, zerophase=True) #Good for 2023-09-02
#    tr_filt.filter('bandpass', freqmin=0.002, freqmax= 0.009, corners=2, zerophase=True) #Good for 2023-09-02
    data   = tr_filt.data
    #The sampling rate is 1 second
    date_rng      = pd.date_range(start= Date_of_Day, freq='%ds'%(sampling_rate), periods = data.size)
    #change the date_rng to numpy
    time          = date_rng.to_numpy()
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    #Envelope of filtered data
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    #Plot the amplitude spectrum
    #ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5, label = 'Spectral %s'%(param))
    #ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5, label = 'Envelope')
    ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5, label = AphaDict[ix])
    #Fill the space between the Figure and the x-xais
    #ax.fill_between(time, data_envelope, color = 'r', alpha = 0.3)
    ax.fill_between(time, data_envelope, color = color, alpha = 0.3)
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=16,                 # Font size for the annotation
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
#    locations = verical_grid()
#    ax.xaxis.set_major_locator(locations)
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
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
#   plt.yticks(fontsize = 11)
    #number of vertical lines for grid


#def read_mseed(FILE):
def read_mseed(Date_of_Day, FILE, xmlFile, comp, bandpass, Pressure = None):
    #get the stream by reading the mseed file
    st      = read(FILE)
    #get all the channls in the stream
    channels = [tr.stats.channel for tr in st]  
    #get the channel you want to extract
    chanel = channels[0].replace(channels[0][-1], comp)
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
    #sorted the stream
    st.sort(keys=["channel"], reverse=True)
    #print(channels,chanel) 
    #extract the trace for the choosing channel
    tr       = st.select(channel=chanel)[0]
    #get the parameters of the trace
    network       = tr.stats.network
    station       = tr.stats.station
    #channel       = tr.stats.channel
    #component     = tr.stats.component
    #print("component :"+ component)
    stime         = str(tr.stats.starttime)
    endtime       = str(tr.stats.endtime)
    #remove mean and trend on the trace
    tr.detrend(type='demean')
    #Title  = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
    Title         = "%s.%s.%s"%(network, station, chanel)
    ##Resampling the data
    tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
    #################################
    #sampling_rate = tr.stats.sampling_rate
    data         = tr.data
    #time          = tr.times()
    #The sampling rate is 1 second
    date_rng     = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
    #change the date_rng to numpy
    time         = date_rng.to_numpy()
    return(time, tr, Title)
    #return(time, tr, Title)


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


def Extract_df_list(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
    #Check if the desired velocity is the Horizontal velocity
    if(param=="velocity_up" or param=="velocity_down"):
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
    elif(param=="turbidity" or param=="avg_BS"):
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



def Plot_Temp2D(ax, fig_object, time, data2D, depths, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    vmind         = np.min(data2D, axis = 1)
    vmaxd         = np.max(data2D, axis = 1)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
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
    im            = ax.imshow(data2D, cmap = color, vmin = min(vmind) , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
                   interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
         xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=16,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax  = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    ########################
    #Add the colorbar on the figure
    fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
    #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[6, 6.5,7 ,7.5,8.5,9])
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
    ax.yaxis.set_major_locator(mticker.MultipleLocator(base=20))
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    ax.set_ylim(float(min(depths)), float(max(depths)))
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])

def Extract_matrix(df, Tlist,  param, nrange=None):
        if(nrange== None):
            r        = df.range.data
        else:
            r        = df.range.data[:nrange]
        if(param=="velocity2D_UP" or param=="velocity2D_DOWN"):
            U        = df.velds.U_mag
            #r        = U.range.data[:nrange]
            Matrix2D = np.concatenate([U.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([U.sel(time = ti)['time'] for ti in Tlist])
            #Free memory
            del U
            return (time, r, Matrix2D)
        elif(param=="veldir2D_UP" or  param=="veldir2D_DOWN"):
            #Get the velocity Direction in 2D
            U        = df.velds.U_dir
            #r        = U.range.data[:nrange]
            Matrix2D = np.concatenate([U.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([U.sel(time = ti)['time'] for ti in Tlist])
            #Free memory
            del U
            return (time, r, Matrix2D) 
        elif(param=="NORTH_VEL"):
            #Get the 2D Northen velocity
            #v =====> to the NORTH
            v2d      = df.velds.v
            #r        = v2d.range.data[:nrange]
            v2d_sel  = np.concatenate([v2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([v2d.sel(time = ti)['time'] for ti in Tlist])
            #get 1D Northern velocity
            v1d      = np.nanmean(v2d_sel.data, axis= 0)
            #Free memory
            del v2d
            #retrun value
            return  time, r, v1d, v2d_sel
        elif(param=="EAST_VEL"):
            #Get the 2D Northen velocity
            #u =====> to the EAST
            #r        = u2d.range.data[:nrange]
            u2d      = df.velds.u
            u2d_sel  = np.concatenate([u2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([u2d.sel(time = ti)['time'] for ti in Tlist])
            #get 1D Eastern velocity
            u1d      = np.nanmean(u2d_sel.data, axis= 0)
            #Free memory
            del u2d
            #return values
            return  time, r, u1d, u2d_sel
        elif(param == "HORIZONTAL_VEL"):
            #Get the 2D Vertical velocity
            #w =====> is the vertical component of the velocity
            #r        = U2d.range.data[:nrange]
            U2d      = df.velds.U_mag
            U2d_sel  = np.concatenate([U2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([U2d.sel(time = ti)['time'] for ti in Tlist])
            #get 1D Eastern velocity
            U1d      = np.nanmean(U2d_sel.data, axis= 0)
            #Free memory
            del U2d
            #return values
            return  time,r, U1d, U2d_sel
        elif(param == "HORIZONTAL_VEL"):
            #Get the 2D Vertical velocity
            #w =====> is the vertical component of the velocity
            #r        = U2d.range.data[:nrange]
            U2d      = df.velds.U_mag
            U2d_sel  = np.concatenate([U2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([U2d.sel(time = ti)['time'] for ti in Tlist])
            #get 1D Eastern velocity
            U1d      = np.nanmean(U2d_sel.data, axis= 0)
            #Free memory
            del U2d
            #return values
            return  time,r, U1d, U2d_sel
        elif(param == "3Comps"):
            #Get the 2D Northen velocity
            #r        = u_2d.range.data[:nrange]
            #u =====> to the EAST
            u_2d      = df.velds.u
            u2d_sel  = np.concatenate([u_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            ########################
            #Get the 2D Northen velocity
            #v =====> to the NORTH
            v_2d     = df.velds.v
            #r        = v_2d.range.data[:nrange]
            v2d_sel  = np.concatenate([v_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            #Get the 2D Vertical velocity
            #w =====> is the vertical component of the velocity
            w_2d     = df.velds.w
            #r        = w_2d.range.data[:nrange]
            w2d_sel  = np.concatenate([w_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            ##########################
            time     = np.concatenate([u_2d.sel(time = ti)['time'] for ti in Tlist])
            #return values
            return   u2d_sel, v2d_sel, w2d_sel
        #elif(param=="backscatter_up"):
        elif(param=="backscatter"):
            alpha = 0.6
            #alpha = 0.1
            time, Matrix2D = BSC_Correct(df, Tlist, alpha, nrange, MEAN=False)
            return (time, r, Matrix2D)
        elif(param=="backscatter_down"):
            alpha = 2.4
            #alpha = 1.0
            time, Matrix2D = BSC_Correct(df, Tlist, alpha, nrange, MEAN=False)
            return (time, r, Matrix2D)



#def Extract_matrix(df, Tlist, param):
#        if(param=="velocity2D"):
#            U        = df.velds.U_mag
#            #r        = U.range
#            r        = U.range.data[:11]
#            Matrix2D = np.concatenate([U.sel(time = ti, range = r) for ti in Tlist], axis =1)
#            #Free memory
#            del U
#        elif(param=="backscatter"):
#            #make the average on all the 4 beams, NB the key word here is amp
#            df_beam_avg = np.mean(df.amp, axis = 0)
#            #Get the range, altitude
#            #r           =  df_beam_avg.range
#            r           =  df_beam_avg.range.data[:11]
#            #Loop over the list of the time
#            Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist], axis =1 )
#            #Free memory
#            del df_beam_avg
#        ####################
#        return Matrix2D


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
        data   = np.mean(DATA2D, axis = 0)
        return (TIMES, data)
    else:
        return(TIMES,DEPTHS, DATA2D)



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
    #Compute the average of the horizontal velocity known as U_mag and it direction known as U_dir
    #df_avg["U_mag"] = df_avg.velds.U_mag
    #df_avg["U_dir"] = df_avg.velds.U_dir
    ##Assign the Turbulence  parameter to df_avg
    #df_avg["TI"]    = df_avg.velds.I

    #return df_avg
    return df_avg




def hist_nbin(data):
    nbins = int(1+3.22 * np.log10(data.size))
    #nbins = np.arange(0, nbins +1, 0.15)
    nbins = np.arange(0, nbins +1, 1)
    return(nbins)

#Define function
def Plot_fig(ax, time, data, param, plot_bigger, label=None,  **PARAMSDICT):
    ##Get the parameters
    fsize = 14
    ylabel, color, ix = PARAMSDICT[param]
    #ax.plot(time, data, kwargs)
    #color='g'
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    #Check the user need the plot to be bigger
    if(plot_bigger):
        ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        #ax.plot(time, data, lw=6.0, linestyle ="-", color = color, alpha =0.5,label = AphaDict[ix])
        ax.plot(time, data, lw=6.0, linestyle ="-", color = color, alpha =0.5)
    else:
        #ax.plot(time, data, lw=0.4, linestyle ="-", color = color, alpha =0.9, label = param.title())
        #ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.9, label = param.title())
        if(label):
            ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.8, label = label)
            #ax.legend(loc="upper right",fontsize='x-large',bbox_to_anchor=(1, 1),frameon=False, shadow=False)
            #Second line on secondary Y-axis
            ax2 = ax.twinx()
            #####################
            ax2.set_ylabel(label,fontsize=16, 
                            #fontweight='bold', 
                            rotation=90)
            #ax.legend(loc="upper left")
            #disable axis
            #ax2.set_yticklabels([])
            ax2.tick_params(axis='y', labelright=False, right=False)
        else:
            ax.plot(time, data, lw=0.7, linestyle ="-", color = color, alpha =0.8)

    #Arrange for Lake Level
    #get the coordinate of the point, where to write the text
    #Write the text on the figure
    if(param == "Lake_Level"):
        #################################
        #ylabel = "Lake Level (m a.s.l)"
        #Set y-axis to meters above sea level
        #ax.get_yaxis().get_major_formatter().set_useOffset(False)
        formatter = ticker.ScalarFormatter() 
        formatter.set_useOffset(False)
        ax.yaxis.set_major_formatter(formatter)
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.2f}'.format(val)))

    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


    #get the coordinate of the point, where to write the text
#    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=   xy_coords,                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.85),                   # Coordinates (relative to axes) for the annotation
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
#    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #Add the legend
    #number of vertical lines for grid
#    locations = verical_grid()
#    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #Label size
    ax.tick_params(axis='both', labelsize=Figsize)
    #ax.legend(loc="upper left")
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    plt.yticks(fontsize = fsize)
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
    #ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])



#def Plot_fig2D(ax, fig_object, time,r, data2D, beam_angle, height_adcp,  blank_dist,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
def Plot_fig2D(ax, fig_object, df, time, data2D,  height_adcp,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    ylabel, color, ix  = PARAMSDICT[param]
    vmind         = np.min(data2D, axis = 1) 
    vmaxd         = np.max(data2D, axis = 1)
    if("velocity" in ylabel):
        vmind = 0.0; vmaxd = 0.4
    #####Test sur L'ADCP orientation
    if(ADCP_DOWN):
        #invert the y-axis
        #ax.invert_yaxis()
        #invert ticks values from positive to negative 
        depths     = abs((height_adcp+ blank_dist) - np.cos(np.radians(beam_angle)) * ranges)         
        y_lims     = [min(depths), height_adcp]
        #y_lims     = [0, max(depths)]
    else:
        depths     = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges      
        y_lims     = [height_adcp, max(depths)]
        #y_lims     = [max(depths), 0]
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time) 
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime' 
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime) 
    end_num       = mdates.date2num(enddatetime) 
    #Set some generic y-limits.
    #y_lims        = [0, max(depths)]
    #y_lims        = [min(depths), max(depths)]
    #Set the extent for 2D plot
    extent        = [start_num , end_num,  y_lims[0], y_lims[1]]
    #Make a 2D plot
    im    = ax.imshow(data2D, cmap = color, vmin = min(vmind) , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
                   #resample=True, extent = extent)
                   interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
    
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position() 
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
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
            ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            ax.yaxis.set_label_coords(-0.08 , 1.2)
        #else:
        #    #fig_object.colorbar(im, cax=cbar_ax, ticks =[0, 90, 180, 270, 350])
        #    ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = 12)
    #elif('CW' not in ylabel):
    else:
        if(ADCP_DOWN):
            #set the color bar on the figure 
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.05, cbar_width, cbar_height])
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.1, 0.2, 0.3, 0.4])
            cbar_ax.yaxis.set_label_position("right")
            #position of the color bar
            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
            #ylabel and the postion of its positions 
            #ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = 12)
            ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
            ax.yaxis.set_label_coords(-0.08 , 1.2)
            #plt.yticks(fontsize = 12)
#        else:
#            #fig_object.colorbar(im, cax=cbar_ax)
#            ax.set_ylabel("H (m)", labelpad = pad, loc="top", fontsize = Figsize)
            #cbar_ax.yaxis.set_label_position("right")
    #Set label
    #Control font size 
    #cbar_ax.tick_params(labelsize=10)
    ##############################################################
    #ax.set_yticklabels(ticks)
    #depths_ticks   = ["%.1f"%(it) for it in depths]
    # depths_ticks   = [float("%.1f"%(it)) for it in depths]
    # ax.set_yticklabels(depths_ticks)
    #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
    if(ADCP_DOWN):
        #invert the y-axis
        #ax.invert_yaxis()
        ##ax.set_ylim(0, float(height_adcp))
#        ax.invert_yaxis()
        #get yticks
        yticks = ax.get_yticks()
        #reverse ticklabels
        ytick_labels = yticks[::-1]
        ticks_new = []
        for ticks in ytick_labels: 
            if (float(ticks) > height_adcp):
                ticks = height_adcp
            ticks_new.append(ticks)
        #ax.set_yticklabels(yticks[::-1])
        #Set the ticks and labels
        ax.set_yticklabels(ticks_new)
        # Optionally, adjust the axis limits to match the ticks
        #ax.set_ylim(min(ticks_new), max(ticks_new))
        ax.set_ylim(max(ticks_new), min(ticks_new))
    else:
        #ax.set_ylim(float(min(depths)), float(max(depths)))
        ax.set_ylim(float(height_adcp), float(max(depths)))
        #invert ticks values from positive to negative 
        ax.annotate(
         AphaDict[ix],     # Text for the annotation
         #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
         xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=16,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
    #Set xlim of x-axis
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
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
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])

##############################################

def plot_spectrogram(ax, fig_object, tr, Date_of_Day, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
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
    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    #ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        xy = xy_coords,                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=16,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position() 
    cbar_width  = 0.007
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.007
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #set the colorbar annotation
    #set the color bar limit on the image
    img.set_clim(-200, -100)  # Set color limits
    #fig_object.colorbar(img, cax = cbar_ax, label = 'Decibel (dB)', format="%+2.f", location='right')
    cbar       = fig_object.colorbar(img, cax = cbar_ax, pad=0.0007, label = 'Power (dB)', format="%+2.f")
    image_data = img.get_array()
    PSD_MEANS  = np.mean(image_data, axis = 1)
    #################################
    #ENER_     = np.sum(np.square(spec))
    ENER_      = np.mean(spec, axis =0)
    #ENER_      = np.square(spec)
    ENRER_db   = 10 * np.log10(ENER_)
    #print(spec.shape)
    #print(len(times))
    #print(len(ENRER_db))
    #exit()
#    print("-----"*50)
#    print(PSD_MEANS)
#    print("-----"*50)
#    print(len(ENRER_db), len(PSD_MEANS))
    file_name  = 'SEMIC_PDS_%s.dat'%(Date_of_Day)
    #Save the array into na file
    #np.savetxt(file_name , PDS_MEANS , delimiter="\n", fmt="%.4f")
    np.savetxt(file_name , ENRER_db , delimiter="\n", fmt="%.4f")
    #np.savetxt(file_name ,abs(ENRER_db) , delimiter="\n", fmt="%.4f")
    #fig_object.colorbar(img, cax = cbar_ax, format="%+2.f dB", pad=0.2)
    #Annotate the colorbar
    #####################
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #disable ticks on  x-axis
    ax.set_xticklabels([])
    return cbar

##############################################
def Func_PPSD(ax, fig_object, Date_of_Day, fmseed, fmxl, comp, Pressure= None, **PARAMSDICT):
    #get the parameters
    #get the value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
    #waveform = read("example.mseed")
    stream    = read(fmseed)
    #Load station metadata (StationXML file)
    inventory = read_inventory(fmxl)
    #get all the channls in the stream
    channels = [tr.stats.channel for tr in stream]  
    #get the channel you want to extract
    chanel   = channels[0].replace(channels[0][-1], comp)
    #Select first trace
    #trace     = stream[0]
    #get the sampling interval
    #sorted the stream
    stream.sort(keys=["channel"], reverse=True)
    #print(channels,chanel)
    #extract the trace for the choosing channel
    trace       = stream.select(channel=chanel)[0]
    ######################################
    #Create PPSD object
    ppsd          = PPSD(trace.stats, inventory)
    #Add entire 24-hour waveform data
    ppsd.add(stream)
    #Original data
    psd_matrix = np.array(ppsd.psd_values)       # shape: (n_times, n_periods)
    #Extract the periods from the ppsd
    periods = np.array(ppsd.period_bin_centers)
    #Extract the time from the ppsd
    times = [t.datetime for t in ppsd.times_processed]
    #Convert the times into the matplotlib
    time_nums = mdates.date2num(times)
    #Transpose: shape (n_periods, n_times)
    Z_psd     = psd_matrix.T
    # --- Calculate extended time range: from 00:00 of day before to 00:00 of day after ---
    #first_time = times[0]
    #creates a new time at midnight (00:00) of the next day (24:00 boundary)
    #to Ensures your plot covers a full 24-hour period
    last_time     = times[-1]
    # 00:00 of day before
    #start_time    = first_time.replace(hour=23, minute=59, second=0, microsecond=0) - timedelta(days=1)
    #start_time_num = mdates.date2num(start_time)
    #extended the last value of the times (last_time) to 24:00:00:00
    next_time     = last_time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(days=1)
    #end_time      = last_time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(days=1)
    #Convert the value to matplot time
    next_time_num = mdates.date2num(next_time)
    #end_time_num  = mdates.date2num(end_time)

    #Adds a new axis (np.newaxis) to make it 2D so it can be horizontally stacked later
    #first_column  = Z_psd[:, 0][:, np.newaxis]    # shape: (n_periods, 1)
    last_column   = Z_psd[:, -1][:, np.newaxis]             # shape: (n_periods, 1)

    #Horizontally appends the last column, repeating the last PSD value flat.
    #This is like saying “assume the last value continues until 24:00.”
    Z_extended          = np.hstack([Z_psd, last_column])          # shape: (n_periods, n_times + 1)
    #Z_extended   = np.hstack([first_column, Z_psd, last_column])  # shape: (n_periods, n_times + 2)

    #Extends the time array to include the extra right edge at midnight.
    #This matches the new column added to Z.
    time_nums_extended  = np.append(time_nums, next_time_num)
    #time_nums_extended  = np.concatenate([[start_time_num], time_nums, [end_time_num]])

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
    #set the limit of y-axis and for the image to compare in the same scale
#    img.set_clim(-200, -100)
#    img.set_clim(-165, -90)
#    ax.set_ylim(min(periods), 1e+2 * 7)
    #ax.set_ylim(1e-1 * 10, 1e+2 * 7)
    ax.set_ylim(1e-2 * 6, 1e+2 * 7)
    #ax.set_ylim(None, 1e+2 * 7)
    #Set axis limits explicitly (to avoid showing extra data outside the range)
    ax.set_xlim(time_nums_extended[0], time_nums_extended[-1])
    #######################
#    ax.set_ylim(4, 200)
#    #ax.axhline(y=100, color='k',linewidth=1.5, linestyle='--', alpha=0.3)
#    ax.axhspan(90, 100*1.2, facecolor="k", alpha = 0.2)
    #############################
    #ax.set_xlim(start_time_num, end_time_num)

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
    #cbar_width  = 0.02
    cbar_width  = 0.007
    cbar_height = bbox.height
    #cbar_x      = bbox.x1 + 0.02
    cbar_x      = bbox.x1 + 0.007
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #check if pressure
    if(Pressure != None):
        cbar = fig_object.colorbar(img, cax=cbar_ax,  pad=0.0007, label = r'PPSD $(\mathrm{ dB \ ref\ Pa^2/Hz})$', format="%+2.f")
    else:
        cbar = fig_object.colorbar(img, cax=cbar_ax,  pad=0.0007, label = r'PPSD $(\mathrm{dB \ ref \ m^2s^{-4}/Hz})$', format="%+2.f")
     #ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        #xy=(0.0089, 0.95),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.92),                   # Coordinates (relative to axes) for the annotation
        xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=16,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
    #####################
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #Label size
    ax.tick_params(axis='both', labelsize=Figsize)
    #ax.set_ylim(1e-2, 1e3)
    #ax.set_ylabel([])
    #disable ticks on  x-axis
    ax.set_xticklabels([])
    return cbar





def plot_seismic_power(ax, tr, Date_of_Day, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    #get the value from the Dictionary
    ylabel, color, ix = PARAMSDICT[param]
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
    # Calculate spectrogram
    #freqs, times, spec = spectrogram(data, fs=sampling_rate, nperseg=256, noverlap=128, scaling='density')
    freqs, times, spec = scipy.signal.spectrogram(data, fs=1.0, nperseg=256, noverlap=128)
    ##########
    #get the coordinate of the point, where to write the text
    #which axes object it should be near Calculate the colorbar position and size
    #free the memory from the data
    del data
    #set the colorbar annotation
    #################################
    ENER_           = np.mean(spec, axis =0)
    #Call the interpolation
    time_new, ENER_ = inter1D(times, ENER_, time)
    #############################
    ENRER_db   = 10 * np.log10(ENER_)
    ##Plot Seismic power energy ###########
    ax.plot(time, ENRER_db, lw=2.0, linestyle ="-", color = color, alpha =0.9)
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    #Set ylim of y-axis
    ax.set_ylim(min(ENRER_db), max(ENRER_db))
    #Set xlim of x-axis
    ax.set_xlim(min(time), max(time))
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
#    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         xy=xy_coords,                   # Coordinates (relative to axes) for the annotation
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
    #ax.grid(visible = True, axis = "both", alpha = 0.7)
    ax.grid(visible = True, axis = "y", alpha = 0.5)
    ax.grid(visible = True, axis = "x", alpha = 0.4)
    #disable ticks on  x-axis
    ax.set_xticklabels([])


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


##############plot histogram ######################
def Plot_hist(ax, time, data, param, plot_bigger,  **PARAMSDICT):
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    nbins             = hist_nbin(data)
    #print(nbins)
    #ax.bar(time, data, width = 0.01, align ='edge')
    #compute the width as function of the data
    #xwd = 1/data.size + (1./data.size) * 0.89
    if(len(set(data)) > 2):
        #xwd = 1/data.size + (1./data.size) * 0.3
        xwd = 1/data.size + (1./data.size) * 0.1
        #ax.bar(time, data, width = (1/data.size) *1.0, align ='edge', color = color, lw= 1000.0, label = param.title())
        ax.bar(time, data, width = xwd , align ='edge', color = color, lw= 100.0, label = param.title())
    else:
        ax.plot(time, data, lw=2.0, linestyle ="-", color = color, alpha =1.0, label = param.title())
    #ax.legend(loc="upper left", prop = {'size': 6})
    #Annotate the figure
    ax.annotate(
        AphaDict[ix],     # Text for the annotation
        xy=   xy_coords,                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.90),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.85),                   # Coordinates (relative to axes) for the annotation
        #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
        xycoords='axes fraction',    # Use axes fraction for positioning
        fontsize=16,                 # Font size for the annotation
        #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
        bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    ###################
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    ###################################
    #Set ylim of x-axis
    ax.set_xlim(tstart, tend)
    #Set label
    ax.set_ylabel(ylabel,  fontsize = Figsize)



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
#with open("configs.yaml") as Fym:
with open("config_OBS_vs_LandBased.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
PARAMSDICT     = Fig_params['PARAMSDICT']
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFiles_OBS      =   Fig_params['mseedFiles_OBS']
mseedFile_Land     =   Fig_params['mseedFile_Land']
#Seimic mseed file
Response_File_Land = Fig_params['Response_File_Land']
Response_File_OBS  = Fig_params['Response_File_OBS']
#Grab the path of the data
ADCP_FILE_NAME_UP = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
#Discharge file name
FILE_DISCH     = Fig_params['FILE_DISCH']
#get the Lake level FILE
FILE_LAKE      = Fig_params['FILE_LAKE']
#Turbidity file name
FILE_TURB      = Fig_params['FILE_TURB']
#Grab the Meteo data file
FILE_METEO     = Fig_params['FILE_METEO']
#Grab the temperature from BRB-SOLO-3
FILE_TEMPE     = Fig_params['FILE_TEMPE']
#############################################
#Check the user want to make a plot bigger?
plot_bigger    = Fig_params['plot_bigger']
#Grab the starttime
STARTDATE      = Fig_params['STARTDATE']
#Grab the ending time
ENDDATE        = Fig_params['ENDDATE']
#Plot Velocity Parameters it True or False
velocity_up    = Fig_params['velocity_up']
velocity_down  = Fig_params['velocity_down']
#Plot the vertical velocity Parameters it True or False
vertical_vel   = Fig_params['vertical_vel']
#Plot Pressure Parameters it True or False
pressure_up    = Fig_params['pressure_up']
pressure_down  = Fig_params['pressure_down']
#get the wind speed
wind_speed     = Fig_params['wind_speed']
##Grab the Temperature to visualize
temperature_up = Fig_params['temperature_up']
temperature_down= Fig_params['temperature_down']
##Grab the Current direction paramter to plot
veldir2D_UP    = Fig_params['veldir2D_UP']
veldir2D_DOWN  = Fig_params['veldir2D_DOWN']
#Grab the temperature 2D
Temperature2D  = Fig_params['Temperature2D']
#density       = Fig_params['density']
#Get the precipation parameter
precipitation  = Fig_params['precipitation']
#Grab the discharge parameter
discharge      = Fig_params['discharge']
############################
Lake_Level     = Fig_params['Lake_Level'] 
#Plot the Turbidity?
turbidity      = Fig_params['turbidity']
#Plot average backscatter
avg_BS         = Fig_params['avg_BS']
##Plot the backscatter of Particules ###
backscatter    = Fig_params['backscatter']
#Tracking surface
surface_tracking = Fig_params['surface_tracking']
#plot the wind direction
wind_direction  =    Fig_params['wind_direction']
#Plot the 2D velocity
velocity2D_UP   = Fig_params['velocity2D_UP']
velocity2D_DOWN = Fig_params['velocity2D_DOWN']
#Plot the seismic Waveform 
waveform        = Fig_params['waveform']
#Components of the seismograms you want to plot
Components       = Fig_params['Components']
pressure_waveform = Fig_params['pressure_waveform']
#plot the seismic power
seismic_power     = Fig_params['seismic_power']
#Plot the seismic Envelope
envelope         = Fig_params['envelope']
#Plot seismic spectrogram
spectrogram    = Fig_params['spectrogram']
#plot the probabilistic density
PPSD_spec      = Fig_params['PPSD_spec']
#Grab the bandpass
bandpass       = Fig_params['bandpass']
bandpass_spec  = Fig_params['bandpass_spec']
#get the nbin
n_bin          = Fig_params['n_bin']
######################################
ADCP_DOWN     = Fig_params['ADCP_DOWN']
ADCP_UP       = Fig_params['ADCP_UP']
########################################
ADCP_DOWN_HEIGTH_FROM_LAKEBED       = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']
ADCP_UP_HEIGTH_FROM_LAKEBED       = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
#This allow to plot the temperature and the pression on the same axis
#Grab the time list
try:
    date_all = pd.date_range(start= STARTDATE, end = ENDDATE)
except:
     print("Checked the date entring %s  %s"%(STARTDATE, ENDDATE))
     exit()
#Create a time list for the started and the
date_list = [str(it).split()[0] for it in date_all]
################################# Provide the files ###################

#mseedFile_OBS ={os.basename(f) : f for f in mseedFiles_OBS}
mseedFile_OBS =dict() 
for f in mseedFiles_OBS:
    basename = os.path.basename(f).split('_')[-1].split('-')[0]
    mseedFile_OBS[basename] = f
###########

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

#index_axis   = [l[-1] for l in PARAMSDICT.values()]
#print(PARAMSDICT)
#print(index_axis)
#remove the wavefrom in the dictionary
waveform_item = PARAMSDICT.pop('waveform')
#Now set the component of the landbased stations
#components list
#Dictionary keys list
List_Keys    = [item for item in PARAMSDICT if 'Atm_pressure' in item or 'wind_speed' in item or 'Lake_Level' in item]
#extend the keys list
List_Keys.extend(Components)
for k in PARAMSDICT:
    if(k not in List_Keys):
        List_Keys.append(k)

#Create a new dictionary
#PARAMSDICT_NEW  = dict()

for ix, item in zip(range(len(List_Keys)), List_Keys):
    if item in PARAMSDICT:
        # Change the axis index
        PARAMSDICT[item][-1] = ix
        #PARAMSDICT_NEW[item] = PARAMSDICT[item]
        PARAMSDICT[item] = PARAMSDICT[item]
    else:
        # Make a new copy of waveform_item and update its index
        #new_waveform = ['µm/s²', 'k', ix]  # Fresh copy each time
        new_waveform = [waveform_item[0], waveform_item[1], ix]  # Fresh copy each time
        PARAMSDICT[item] = new_waveform
        #PARAMSDICT_NEW[item] = new_waveform
        PARAMSDICT[item] = new_waveform

#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs        = len(PARAMSDICT)


#Read the ADCP file
try:
    if(ADCP_UP):
        #read ADCP looking upward
        #df      = dlfn.read(ADCP_FILE_NAME_UP)
        EXT     = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1]
        #########################################
    elif(ADCP_DOWN):
        #df      = dlfn.read(ADCP_FILE_NAME_DOWN)
        EXT     = os.path.basename(ADCP_FILE_NAME_DOWN).split("_")[1]
    else:
        #df      = dlfn.read(ADCP_FILE_NAME_UP)
        EXT     = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1]
############################################
except:
    print("Provide the two  ADCP-Files that need to be plot")
    exit()

#***Just for pressure**********#
df_up          = dlfn.read(ADCP_FILE_NAME_UP)
df_avg_up      = Average(df_up, n_bin)
beam_angle_up  = df_avg_up.beam_angle
blank_dist_up  = df_avg_up.blank_dist
#####*****************#######

#Read the ADCP file
#df      = dlfn.read(ADCP_FILE_NAME) 
df      = dlfn.read(ADCP_FILE_NAME_DOWN)
######Perform the averaging#
df_avg = Average(df, n_bin)
#get the beam angle
beam_angle   = df_avg.beam_angle
#get the blink distance
blank_dist   = df_avg.blank_dist
############################################
ranges       = df_avg.range.data
#Create the figure




#Read the Turbidity  file
d_echo = pd.read_csv(FILE_TURB, delimiter='\t')
#Add the of DateTime
d_echo = Add_Column(d_echo)
#Read the Turbidity  file
#d_disch= pd.read_fwf(FILE_DISCH, delimiter='\t')
##Add the of DateTime
#d_disch= Add_Column(d_disch)

#load the Lake Level Data
dfLake = pd.read_csv(FILE_LAKE, skiprows=0, sep =",",encoding="latin1")
#Load the Discahge file
d_disch= pd.read_csv(FILE_DISCH, skiprows=0, sep =",",encoding="latin1") 
#print(d_disch) 
#read the Meteo database
d_mteo = pd.read_csv(FILE_METEO, delimiter ='\t')
#Add the of DateTime
d_mteo = Add_Column(d_mteo)

#Grab 2D Temperature
#get the temperature data from RBR-SOLO-3
# Opening JSON file
ftmp   = open(FILE_TEMPE)
# returns JSON object as
# a dictionary
data_temp = json.load(ftmp)
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])


#Create the figure
#fig, axs= plt.subplots(nfigs, 1, sharex = False, figsize=(12,10))
fig, axs= plt.subplots(nfigs, 2, sharex = False, figsize= fig_size)
#set Dictionnary for components OBS vs Landbased

#OBS_comps = {'Z': '1', 'N':'2', 'E': '3'}
OBS_comps = {'Z': 'HH1', 'N':'HH2', 'E': 'HH3'}

#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#fig, axs  = plt.subplots(6, 1, sharex = False, figsize = fig_size)
#fig, axs  = plt.subplots(nfigs+2, 1, sharex = False, figsize = fig_size)
#### Plot the envolop ##############
#loop over the axis

if(envelope):
    #set the parameter
    param            = "envelope"
    #get the seismogram from mseedFile
    time_obs, tr_obs, Title_obs     = read_mseed(STARTDATE, mseedFile_OBS, Response_File_OBS, bandpass_spec)
    time_Land, tr_Land, Title_Land  = read_mseed(STARTDATE, mseedFile_Land, Response_File_Land, bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #set the title
    #axs[ix].set_title(Title, loc='center', pad=None)
    axs[x_i, y_j].set_title(Title, loc='center', pad=None)
    #Plot the figure by calling the plotting
    #PlotEnvelop(axs[ix], STARTDATE, param, tr, linewidth=1.2,  **PARAMSDICT)
    PlotEnvelop(axs[x_i, y_j], STARTDATE, param, tr, linewidth=1.2,  **PARAMSDICT)
#plot seismic waveform
if(waveform):
    for comp in Components:
        #set the parameter
        param      = comp
        #get the seismogram from mseedFile
        time_Land, tr_Land, Title_Land  = read_mseed(STARTDATE, mseedFile_Land, Response_File_Land, comp, bandpass_spec)
        #get the OBS-component
        comp_obs                    = OBS_comps[comp]
        #get the OBS data
        #print(comp, comp_obs)
        time_obs, tr_obs, Title_obs = read_mseed(STARTDATE, mseedFile_OBS[comp_obs], Response_File_OBS, comp_obs[-1], bandpass_spec)
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix               = PARAMSDICT[param]
        #label        = Title
        for y_j in range(2):
            #Plot the landbased data
            #data  =  tr_filt.data * 1e+6
            ############################################################################
            if(y_j==0):
                Plot_fig(axs[ix, y_j], time_Land, tr_Land.data * 1e+6, param, False,  **PARAMSDICT)
            #tr_filt.filter('highpass', freq=0.2, corners=2, zerophase=True)
            #data  =  tr.data * 1e+6
            ############################################################################
            if(y_j==1):
                Plot_fig(axs[ix, y_j], time_obs, tr_obs.data * 1e+6, param, False,  **PARAMSDICT)
                #remove the y-label
                axs[ix, y_j].set_ylabel(None)
    #Free memory for the time and data
    del time_Land, tr_Land, time_obs, tr_obs
#######################################################
if(pressure_waveform):
    #set the parameter
    param           = "pressure_waveform"
    #get the seismogram from mseedFile
    time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass, Pressure = True)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix   = PARAMSDICT[param]
    #Set the Title
    label           = Title
    axs[ix].set_title(label, loc='center', pad=None)
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, tr.data, param, False, **PARAMSDICT)

if(seismic_power):
        #set the parameter
        param            = "seismic_power"
        #get the seismogram from mseedFile
        #Plot the seismic energy
        #plot_seismic_power(axs[ix], tr, STARTDATE, **PARAMSDICT)
        if 'Z' in Components:
            comp = [it for it in Components if(it=="Z")][0]
        elif('N' in Components):
            comp = [it for it in Components if(it=="N")][0]
        else:
            comp = [it for it in Components if(it=="E")][0]
        #set the parameter
        #param      = comp
        #get the seismogram from mseedFile
        time_Land, tr_Land, Title_Land  = read_mseed(STARTDATE, mseedFile_Land, Response_File_Land, comp, bandpass_spec)
        #get the OBS-component
        comp_obs                        = OBS_comps[comp]
        #get the OBS data
        #print(comp, comp_obs)
        time_obs, tr_obs, Title_obs = read_mseed(STARTDATE, mseedFile_OBS[comp_obs], Response_File_OBS, comp_obs[-1], bandpass_spec)
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix               = PARAMSDICT[param]
        #label        = Title
        for y_j in range(2):
            #Plot the landbased data
            #data  =  tr_filt.data * 1e+6
            ############################################################################
            if(y_j==0):
                plot_seismic_power(axs[ix, y_j], tr_Land, STARTDATE,  **PARAMSDICT)
                #plot_seismic_power(axs[ix], tr, STARTDATE, **PARAMSDICT)
            #tr_filt.filter('highpass', freq=0.2, corners=2, zerophase=True)
            #data  =  tr.data * 1e+6
            ############################################################################
            if(y_j==1):
                #Plot_fig(axs[ix, y_j], time_obs, tr_obs.data * 1e+6, param, False,  **PARAMSDICT)
                plot_seismic_power(axs[ix, y_j], tr_obs, STARTDATE,  **PARAMSDICT)
                #remove the y-label
                axs[ix, y_j].set_ylabel(None)
        #Free memory for the time and data
        del time_Land, tr_Land, time_obs, tr_obs
if(spectrogram):
    #set the parameter
    param            = "spectrogram"
    #get the seismogram from mseedFile
    #############################################################
    if 'Z' in Components:
        comp = [it for it in Components if(it=="Z")][0]
    elif('N' in Components):
        comp = [it for it in Components if(it=="N")][0]
    else:
        comp = [it for it in Components if(it=="E")][0]
    #set the parameter
    #param      = comp
    #get the seismogram from mseedFile
    time_Land, tr_Land, Title_Land  = read_mseed(STARTDATE, mseedFile_Land, Response_File_Land, comp, bandpass_spec)
    #get the OBS-component
    comp_obs                        = OBS_comps[comp]
    #get the OBS data
    #print(comp, comp_obs)
    time_obs, tr_obs, Title_obs = read_mseed(STARTDATE, mseedFile_OBS[comp_obs], Response_File_OBS, comp_obs[-1], bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix               = PARAMSDICT[param]
    #label        = Title
    for y_j in range(2):
        #Plot the landbased data
        ############################################################################
        if(y_j==0):
            ##Plot the figure by calling the plotting function of 2D Matrix
            cbar = plot_spectrogram(axs[ix, y_j], fig, tr_Land, STARTDATE, **PARAMSDICT)
            cbar.remove()
        ############################################################################
        if(y_j==1):
            ##Plot the figure by calling the plotting function of 2D Matrix
            cbar =plot_spectrogram(axs[ix, y_j], fig, tr_obs, STARTDATE, **PARAMSDICT)
            #remove the y-label
            axs[ix, y_j].set_ylabel(None)
            axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)


    #Free memory for the time and data
    del time_Land, tr_Land, time_obs, tr_obs

# PPSD_spec    
if(PPSD_spec):
    #set the parameter
    param            = "PPSD_spec"
    #get the seismogram from mseedFile
    #############################################################
    if 'Z' in Components:
        comp = [it for it in Components if(it=="Z")][0]
    elif('N' in Components):
        comp = [it for it in Components if(it=="N")][0]
    else:
        comp = [it for it in Components if(it=="E")][0]
    #set the parameter
    #param      = comp
    #get the seismogram from mseedFile
    #time_Land, tr_Land, Title_Land  = read_mseed(STARTDATE, mseedFile_Land, Response_File_Land, comp, bandpass_spec)
    #################
    #get the OBS-component
    comp_obs                        = OBS_comps[comp]
    #get the OBS data
    #print(comp, comp_obs)
    #time_obs, tr_obs, Title_obs = read_mseed(STARTDATE, mseedFile_OBS[comp_obs], Response_File_OBS, comp_obs[-1], bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix               = PARAMSDICT[param]
    #label        = Title
    for y_j in range(2):
        #Plot the landbased data
        ############################################################################
        if(y_j==0):
            ##Plot the figure by calling the plotting function of 2D Matrix
            #cbar = plot_spectrogram(axs[ix, y_j], fig, tr_Land, STARTDATE, **PARAMSDICT)
            cbar = Func_PPSD(axs[ix, y_j], fig, STARTDATE, mseedFile_Land, Response_File_Land, comp, **PARAMSDICT)
            cbar.remove()
        ############################################################################
        if(y_j==1):
            ##Plot the figure by calling the plotting function of 2D Matrix
            #cbar =plot_spectrogram(axs[ix, y_j], fig, tr_obs, STARTDATE, **PARAMSDICT)
            cbar=Func_PPSD(axs[ix, y_j], fig, STARTDATE, mseedFile_OBS[comp_obs], Response_File_OBS, comp_obs[-1], **PARAMSDICT)
            #remove the y-label
            axs[ix, y_j].set_ylabel(None)
            #axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
            #keep log gridlines for visual reference
            #axs[ix, y_j].yaxis.grid(True, which='both')
            axs[ix, y_j].set_yticklabels([])
#plot the wind velocity
if(wind_speed):
    #get the Dischage data for the corresponding time
    param         = "wind_speed"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    #label        = Title
    for y_j in range(2):
       #Plot the landbased data
       ############################################################################
       if(y_j==0):
           Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
           ############################################################################
       if(y_j==1):
           Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
           #remove the y-label
           axs[ix, y_j].set_ylabel(None)
           axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
    #Free memory for the time and data
    del time, data
##Plot the Discharge
if(discharge):
        #set the parameter
        param          = "discharge"
        #get the Dischage data for the corresponding time
        #time, data    = Extract_database(d_disch,  date_list, 'Discharge')
        time, data     = ExtractLakeLevel(d_disch,  date_list)
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix  = PARAMSDICT[param]
        #Plot the figure by calling the plotting function, plot_twinx
        #Define the second axis to plot the Lake Level
        for y_j in range(2):
           #Plot the landbased data
           ############################################################################
           if(y_j==0):
               Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
               ############################################################################
           if(y_j==1):
               Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
               #remove the y-label
               axs[ix, y_j].set_ylabel(None)
               axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
##Plot the Lake Level
if(Lake_Level):
        #set the parameter
        param          = "Lake_Level"
        #get the Dischage data for the corresponding time
        #time, data    = Extract_database(d_disch,  date_list, 'Discharge')
        time, data     = ExtractLakeLevel(dfLake,  date_list)
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix  = PARAMSDICT[param]
        #Plot the figure by calling the plotting function, plot_twinx
        for y_j in range(2):
           #Plot the landbased data
           ############################################################################
           if(y_j==0):
               Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
               ############################################################################
           if(y_j==1):
               Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
               #remove the y-label
               axs[ix, y_j].set_ylabel(None)
               axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
##Doubleched if we need to plot the Velocity
if(velocity_up):
   #set the parameter
   param      = "velocity_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
   #Plot the figure by calling the plotting function
if(velocity_down):
   #set the parameter
   param      = "velocity_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg , date_list, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
   #Plot the figure by calling the plotting function
if(vertical_vel):
       #set the parameter
       param      = "vertical_vel"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg , date_list, param)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix = PARAMSDICT[param]
       #Plot the figure by calling the plotting function, plot_twinx
       for y_j in range(2):
          #Plot the landbased data
          ############################################################################
          if(y_j==0):
              Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
              ############################################################################
          if(y_j==1):
              Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
              #remove the y-label
              axs[ix, y_j].set_ylabel(None)
              axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
##Doubleched if we need to plot the Temperature
if(temperature_up):
   #set the parameter
   param        = "temperature_up"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_up , date_list, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix= PARAMSDICT[param]
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)

if(temperature_down):
   #set the parameter
   param        = "temperature_down"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg , date_list, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix= PARAMSDICT[param]
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
##Plot the Pressure
if(pressure_up):
   #set the pad_echr
   param         = "pressure"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg_up , date_list, param)
   #time, data    = Extract_df_list(df_avg , date_list, param)
   param         = param+"_up"
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)

if(pressure_down):
   #set the pad_echr
   param         = "pressure"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg , date_list, param)
   #time, data    = Extract_df_list(df_avg , date_list, param)
   param         = param+"_down"
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_fig(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)
###################################################
if(Temperature2D):
   param            = "Temperature2D"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, depths, data2D       = Extract_RBR_SOLO_TEMP(data_temp,  date_list, AVERAGE=False)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix            = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot the figure by calling the plotting function of 2D Matrix
   #Plot the figure by calling the plotting function, plot_twinx
   for y_j in range(2):
      #Plot the landbased data
      ############################################################################
      if(y_j==0):
          Plot_Temp2D(axs[ix, y_j], fig, time, data2D, depths,  **PARAMSDICT)
          ############################################################################
      if(y_j==1):
          Plot_Temp2D(axs[ix, y_j], fig, time, data2D, depths,  **PARAMSDICT)
          #remove the y-label
          axs[ix, y_j].set_ylabel(None)
          axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)


##Plot the backscatter
if(backscatter):
    #set the parameter
    #param      = "backscatter"
    param      = "backscatter_down"
    #get the Turbidity data for the corresponding time
    #data2D     =  Extract_matrix(df_avg, date_list, param)
    time, r, data2D = Extract_matrix(df_avg, date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
#Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist,ADCP_DOWN=True, rframe= "bottom", **PARAMSDICT)
    Plot_fig2D(axs[ix], fig, df_avg, time, data2D,ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="bottom" , **PARAMSDICT)
##Plot the 2D velocity
if(velocity2D_UP):
    #set the parameter
    param      = "velocity2D_UP"
    #get the Turbidity data for the corresponding time
    #data2D     =  Extract_matrix(df_avg, date_list, param)
    time, r, data2D = Extract_matrix(df_avg_up, date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig,time,r,data2D,beam_angle_up,ADCP_UP_HEIGTH_FROM_LAKEBED,blank_dist_up,ADCP_DOWN=False, rframe= "top", **PARAMSDICT)
    Plot_fig2D(axs[ix], fig, df_avg_up, time, data2D, ADCP_UP_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="top" , **PARAMSDICT)
##Plot the 2D velocity
if(velocity2D_DOWN):
    #set the parameter
    param      = "velocity2D_DOWN"
    #get the Turbidity data for the corresponding time
    #data2D     =  Extract_matrix(df_avg, date_list, param)
    time, r, data2D = Extract_matrix(df_avg, date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist,ADCP_DOWN=True, rframe= "top", **PARAMSDICT)
    Plot_fig2D(axs[ix], fig, df_avg, time, data2D, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=True, rframe="top" , **PARAMSDICT)
##############################
if(veldir2D_UP):
   #set the parameter
   param            = "veldir2D_UP"
   #get the Turbidity data for the corresponding time
   #time, data2D     =  Extract_matrix(df_avg_up, date_list_UTC, "veldir2D", num_range=30)
   #num_range=30
   time, r, data2D  =  Extract_matrix(df_avg_up, date_list, "veldir2D_UP")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix    = PARAMSDICT[param]
   ##############################
   #Plot the figure by calling the plotting function of 2D Matrix
   #Plot_fig2D(axs[ix],fig,time,r,data2D, beam_angle_up, ADCP_UP_HEIGTH_FROM_LAKEBED, blank_dist_up,ADCP_DOWN=False, rframe= "top", **PARAMSDICT)
   Plot_fig2D(axs[ix], fig, df_avg_up, time, data2D, ADCP_UP_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=False, rframe="top" , **PARAMSDICT)
   ###########
if(veldir2D_DOWN):
   #set the parameter
   param              = "veldir2D_DOWN"
   time, r, data2D    =  Extract_matrix(df_avg, date_list, "veldir2D_DOWN")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix      = PARAMSDICT[param]
   #Plot the figure by calling the plotting function of 2D Matrix
   #Plot_fig2D(axs[ix],fig,time,r,data2D,beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED,blank_dist,ADCP_DOWN=True, rframe= "top", **PARAMSDICT)
   Plot_fig2D(axs[ix], fig, df_avg, time, data2D, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  ADCP_DOWN=True, rframe="top" , **PARAMSDICT)

#Plot the 1D backscatter
if(avg_BS):
    #set the parameter
    param      = "avg_BS"
    #get the Turbidity data for the corresponding time
    #time, data = Extract_database(d_echo,  date_list, 'Echo')
    time, data = Extract_df_list(df_avg,  date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
##Plot the Turbidity
if(turbidity):
    #set the parameter
    param      = "turbidity"
    #get the Turbidity data for the corresponding time
    #time, data = Extract_database(d_echo,  date_list, 'Echo')
    time, data = Extract_df_list(df_avg,  date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

#Plot the precipitation
if(precipitation):
    #set the parameter
    param         = "precipitation"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    # convert hPa to bar;  1hPa =   0.001 bar
    #Plot_hist(axs[ix], time, data, param, plot_bigger,  **PARAMSDICT)
    #Define the second axis to plot the Lake Level
    for y_j in range(2):
       #Plot the landbased data
       ############################################################################
       if(y_j==0):
           Plot_hist(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
           ############################################################################
       if(y_j==1):
           Plot_hist(axs[ix, y_j], time, data, param, plot_bigger,  **PARAMSDICT)
           #remove the y-label
           axs[ix, y_j].set_ylabel(None)
           axs[ix, y_j].tick_params(axis='y', labelleft=False)  # Remove y-axis tick labels (graduations)

#Set the axis-label for the last axis
#format='%y %m %d %H %M %S'
#axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%m-%d:%H"))
#axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%H-%M:%S"))

#axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%H:%M"))

#basename_ = os.path.basename(mseedFile) 
#
#basename = basename_.split(".mseed")[0]

basename = 'Fig'


#### Write on the x-axis
if(STARTDATE != ENDDATE):
#    axs[-1].set_xlabel('Time from %s to %s'%(STARTDATE, ENDDATE), fontsize = 13)
    #figname = "ADCP_WAVEFORM_PLOTS_from_%s_to_%s.png"%(STARTDATE, ENDDATE)
    figname = "ADCP_%s_from_%s_to_%s.png"%(basename, STARTDATE, ENDDATE)
if(STARTDATE == ENDDATE):
    #axs.flatten()[-1].xaxis.set_major_locator(HourLocator())
    axs[-1,0].xaxis.set_major_formatter(DateFormatter('%H:%M'))
    axs[-1,0].set_xlabel('Time (hour:minute) on %s'%(STARTDATE), fontsize = 6)
    ################################################
    #axs.flatten()[-2].xaxis.set_major_locator(HourLocator())
    axs[-1, 1].xaxis.set_major_formatter(DateFormatter('%H:%M'))
    axs[-1, 1].set_xlabel('Time (hour:minute) on %s'%(STARTDATE), fontsize = 6)
    figname   = "ADCP_%s_of_%s_TEST.png"%(basename, STARTDATE)
    #figname   = "ADCP_%s_of_%s_TEST_HH.png"%(basename, STARTDATE)
    #figname   = "ADCP_%s_of_%s.pdf"%(basename, STARTDATE)
# Set the font size of yticks
#plt.yticks(fontsize=13)
## Set the font size of xticks
#plt.xticks(fontsize=14)
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
