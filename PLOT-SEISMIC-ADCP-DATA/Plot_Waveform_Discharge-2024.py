#!/usr/bin/env python

#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
import obspy
from obspy import read, read_inventory
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
##############################################
import matplotlib.dates as dt
import xarray as xr
import pandas as pd
from mhkit import dolfyn
from mhkit import dolfyn as dlfn
from datetime import datetime
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




AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}
#size for anotation
f_size = 13
pad = 10
############################################
def InterSpline(t_in, d_in):
    #defined time step
    #tstep = 1.0
    tstep = 0.5
    spl   = InterpolatedUnivariateSpline(t_in, d_in)
    #Interpolation
    t_new = np.arange(min(t_in), max(t_in), tstep)
    d_new = spl(t_new)
    #define the 2D dimension array to store the data
    dim   = (t_new.shape[0], 2) 
    #Define the 2D array
    data_out = np.zeros(dim)
    #Fill the 2D array
    data_out[:,0] = t_new
    data_out[:,1] = d_new
    return data_out

def inter1D(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
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

def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations

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
    #Envelope of filtered data
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    #Plot the amplitude spectrum
    ax.plot(time, data_envelope, lw = linewidth, linestyle ="-", color = color, alpha =0.5, label = 'Spectral %s'%(param))
    #Fill the space between the Figure and the x-xais
    #ax.fill_between(time, data_envelope, color = 'r', alpha = 0.3)
    ax.fill_between(time, data_envelope, color = color, alpha = 0.3)
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    #ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)

    #Remove the frame on the plot
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ##################################
    #if(float(ymax) < 0.01):
#        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    ax.legend(loc="upper left")
    #Set label
    #disable ticks on  x-axis
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    #Remove the ticks marks (small vertical lines) on the x-axis
    ax.tick_params(axis="both", length = 0, color= "white", width = 0)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
#   plt.yticks(fontsize = 11)
    #number of vertical lines for grid


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

def Extract_df_list(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
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
    elif(param=="turbidity"):
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
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
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
    #fig_object.colorbar(im, cax=cbar_ax)
    fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
    #################################
    #Set label
    ax.set_ylabel('H (m)', labelpad = pad, fontsize = 11)
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
    ax.set_yticklabels(ticks)
    #ax.yaxis.set_major_locator(mticker.MultipleLocator(base=5))
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
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
        if(param=="velocity2D"):
            U        = df.velds.U_mag
            #r        = U.range.data[:nrange]
            Matrix2D = np.concatenate([U.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([U.sel(time = ti)['time'] for ti in Tlist])
            #Free memory
            del U
            return (time, r, Matrix2D)
        elif(param=="veldir2D"):
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
            return (time, Matrix2D)
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
def Plot_fig(ax, time, data, param, plot_bigger,  **PARAMSDICT):
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    #ax.plot(time, data, kwargs)
    #Check the user need the plot to be bigger
    if(plot_bigger):
        ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        ax.plot(time, data, lw=6.0, linestyle ="-", color = color, alpha =0.5,label = param.title())
    else:
        #ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0, label = param.title())
        ax.plot(time, data, lw=0.4, linestyle ="-", color = color, alpha =0.9, label = param.title())

    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #Add the legend
#    ax.legend(loc="upper left")
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    plt.yticks(fontsize = 11)
    #Set ylim of x-axis
    ax.set_xlim(min(time), max(time))
    #disable axis
    ax.set_xticklabels([])



def Plot_fig2D(ax, fig_object, time,r, data2D, beam_angle, height_adcp,  blank_dist,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix  = PARAMSDICT[param]
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
    #y_lims        = [0, nrange]
    y_lims        = [0, len(r)]
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
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
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
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
            cbar_ax.yaxis.set_label_position("right")
            #shift the label by 3.5 along x-axis and by 1.0 along the y-axis
            #cbar_ax.yaxis.set_label_coords(3.5, 1.1)
            #move the label by -0.08 on x-axis and by 1.3 on y-axis
            #ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = 12)
            ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = 12)
            #ax.yaxis.set_label_coords(-0.08 , 1.3)
        else:
            fig_object.colorbar(im, cax=cbar_ax, ticks =[0, 90, 180, 270, 350])
    #elif('CW' not in ylabel):
    else:
        if(ADCP_DOWN):
            #set the color bar on the figure 
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            cbar_ax.yaxis.set_label_position("right")
            #position of the color bar
            #cbar_ax.yaxis.set_label_coords(3.5, 1.1)
            #ylabel and the postion of its positions 
            #ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = 12)
            ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = 12)
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
            #ax.yaxis.set_label_coords(-0.08 , 1.3)
            #plt.yticks(fontsize = 12)
        else:
            fig_object.colorbar(im, cax=cbar_ax)
            ax.set_ylabel("H (m)", labelpad = pad, loc="top", fontsize = 12)
            cbar_ax.yaxis.set_label_position("right")
    #Set label
    ##############################################################
    if(ADCP_DOWN):
        #invert the y-axis
        ax.invert_yaxis()
        #invert ticks values from positive to negative 
        depths     = (height_adcp+ blank_dist) - np.cos(beam_angle) * ticks         
    else:
        depths     = (height_adcp + blank_dist) + np.cos(beam_angle) * ticks         
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    #Set the inverted ticks label
    #ax.set_yticklabels(ticks)
    depths_ticks   = ["%.1f"%(it) for it in depths]
    ax.set_yticklabels(depths_ticks)
    #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
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
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position() 
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax     = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    #set the colorbar annotation
    #fig_object.colorbar(img, cax = cbar_ax, label = 'Decibel (dB)', format="%+2.f", location='right',ticks=levels, spacing='proportional')
    fig_object.colorbar(img, cax = cbar_ax, label = 'Decibel (dB)', format="%+2.f", location='right')
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
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
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
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    #Set ylim of y-axis
    ax.set_ylim(min(ENRER_db), max(ENRER_db))
    #Set xlim of x-axis
    ax.set_xlim(min(time), max(time))
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    #ax.grid(visible = True, axis = "both", alpha = 0.7)
    ax.grid(visible = True, axis = "y", alpha = 0.5)
    ax.grid(visible = True, axis = "x", alpha = 0.4)
    #disable ticks on  x-axis
    ax.set_xticklabels([])



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
    ax.legend(loc="upper left", prop = {'size': 6})
    #Set label
    ax.set_ylabel(ylabel,  fontsize = 11)



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
with open("configs.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
PARAMSDICT     = Fig_params['PARAMSDICT']
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFile      =   Fig_params['mseedFile']
#Seimic mseed file
Response_File  = Fig_params['Response_File']
#Grab the path of the data
ADCP_FILE_NAME_UP = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
#Discharge file name
FILE_DISCH     = Fig_params['FILE_DISCH']
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
velocity       = Fig_params['velocity']
#Plot the vertical velocity Parameters it True or False
vertical_vel   = Fig_params['vertical_vel']
#Plot Pressure Parameters it True or False
pressure       = Fig_params['pressure']
##Grab the Velocity component to visualize
temperature    = Fig_params['temperature']
##Grab the Current direction paramter to plot
veldir2D       = Fig_params['veldir2D']
#Grab the temperature 2D
Temperature2D  = Fig_params['Temperature2D']
#density       = Fig_params['density']
#Get the precipation parameter
precipitation  = Fig_params['precipitation']
#Grab the discharge parameter
discharge      = Fig_params['discharge']
#Plot the Turbidity?
turbidity      = Fig_params['turbidity']
##Plot the backscatter of Particules ###
backscatter = Fig_params['backscatter']
#Tracking surface
surface_tracking = Fig_params['surface_tracking']
#plot the wind direction
wind_direction  =    Fig_params['wind_direction']
#Plot the 2D velocity
velocity2D = Fig_params['velocity2D']
#Plot the seismic Waveform 
waveform          = Fig_params['waveform']
pressure_waveform = Fig_params['pressure_waveform']
#plot the seismic power
seismic_power     = Fig_params['seismic_power']
#Plot the seismic Envelope
envelope         = Fig_params['envelope']
#Plot seismic spectrogram
spectrogram    = Fig_params['spectrogram']
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

#Read the ADCP file
try:
    if(ADCP_UP):
        #read ADCP looking downward
        #read ADCP looking upward
        df      = dlfn.read(ADCP_FILE_NAME_UP)
        EXT     = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1]
        #########################################
    elif(ADCP_DOWN):
        df      = dlfn.read(ADCP_FILE_NAME_DOWN)
        EXT     = os.path.basename(ADCP_FILE_NAME_DOWN).split("_")[1]
    else:
        df      = dlfn.read(ADCP_FILE_NAME_UP)
        EXT     = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1]
############################################
except:
    print("Provide the two  ADCP-Files that need to be plot")
    exit()

#***Just for pressure**********#
df_up    = dlfn.read(ADCP_FILE_NAME_UP)
df_up_av = Average(df_up, n_bin)
#beam_up  = df_up_av.beam_angle
#blank_up = df_up_av.blank_dist
#####*****************#######

#Read the ADCP file
#df      = dlfn.read(ADCP_FILE_NAME) 
######Perform the averaging#
df_avg = Average(df, n_bin)
#Create the figure

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
d_disch= pd.read_fwf(FILE_DISCH, delimiter='\t')
#Add the of DateTime
d_disch= Add_Column(d_disch)
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

fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)


#### Plot the envolop ##############
if(envelope):
    #set the parameter
    param            = "envelope"
    #get the seismogram from mseedFile
    time, tr, Title  = read_mseed(STARTDATE, mseedFile, Response_File, bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the figure by calling the plotting
    PlotEnvelop(axs[ix], STARTDATE, param, tr, linewidth=1.2,  **PARAMSDICT)
#plot seismic waveform
if(waveform):
    #set the parameter
    param      = "waveform"
    #get the seismogram from mseedFile
    time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix  = PARAMSDICT[param]
    #Set the Title
    #label        = '%s   %s             BP: %s s'%(basename_new[0], basename_new[1], bandpass)
    #label        = '%s       %s    '%(basename_new[0], basename_new[1])
    #label        = '%s '%(bandpass)
    label        = Title
    if(envelope==False):
        axs[ix].set_title(label, loc='center', pad=None)
    #Plot the figure by calling the plotting function, plot_twinx
    data  =  tr.data * 1e+6
    #time = time.astype("datetime64[ns]")
    Plot_fig(axs[ix], time, data, param, False, **PARAMSDICT)
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
    Plot_fig(axs[ix], time, data, param, False, **PARAMSDICT)

if(seismic_power):
    #set the parameter
    param            = "seismic_power"
    #get the seismogram from mseedFile
    time, tr, Title  = read_mseed(STARTDATE, mseedFile, Response_File, bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the seismic energy
    plot_seismic_power(axs[ix], tr, STARTDATE, **PARAMSDICT)
if(spectrogram):
    #set the parameter
    param            = "spectrogram"
    #get the seismogram from mseedFile
    if(pressure_waveform):
        time, tr, Title  = read_mseed(STARTDATE, mseedFile, Response_File, bandpass_spec, Pressure=True)
    else:
        time, tr, Title  = read_mseed(STARTDATE, mseedFile, Response_File, bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix    = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    plot_spectrogram(axs[ix], fig, tr, STARTDATE, **PARAMSDICT)

##Plot the Discharge
if(discharge):
    #set the parameter
    param         = "discharge"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_disch,  date_list, 'Discharge')
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

##Doubleched if we need to plot the Velocity
if(velocity):
   #set the parameter
   param      = "velocity"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg , date_list, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function
if(vertical_vel):
       #set the parameter
       param      = "vertical_vel"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg , date_list, param)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix = PARAMSDICT[param]
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
       #if(velocity):
       #     plot_twinx(axs[PARAMSDICT['velocity'][-1]] , time, data, param, **PARAMSDICT)
       #     #plot_twinx(axs[ix], time, data, param, **PARAMSDICT)
       #else:
       #     Plot_fig(axs[ix], time, data, param, **PARAMSDICT)
##Doubleched if we need to plot the Temperature
if(temperature):
   #set the parameter
   param        = "temperature"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg , date_list, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix= PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

##Plot the Pressure
if(pressure):
   #set the pad_echr
   param         = "pressure"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_up_av , date_list, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #if(temperature):
   #    plot_twinx(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #else:
   #    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(Temperature2D):
   param            = "Temperature2D"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, depths, data2D       = Extract_RBR_SOLO_TEMP(data_temp,  date_list, AVERAGE=False)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix            = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_Temp2D(axs[ix],fig, time, data2D, depths,  **PARAMSDICT)
   #Plot_Temp2D(axs,fig, time, data2D, depths, **PARAMSDICT)


##Plot the backscatter
if(backscatter):
    #set the parameter
    param      = "backscatter"
    #get the Turbidity data for the corresponding time
    #data2D     =  Extract_matrix(df_avg, date_list, param)
    time, r, data2D = Extract_matrix(df_avg, date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig, data2D, **PARAMSDICT)
    Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist,ADCP_DOWN=True, rframe= "bottom", **PARAMSDICT)

##Plot the 2D velocity
if(velocity2D):
    #set the parameter
    param      = "velocity2D"
    #get the Turbidity data for the corresponding time
    #data2D     =  Extract_matrix(df_avg, date_list, param)
    time, r, data2D = Extract_matrix(df_avg, date_list, param)
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig, data2D, **PARAMSDICT)
    Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist,ADCP_DOWN=True, rframe= "top", **PARAMSDICT)
##############################
if(veldir2D):
   #set the parameter
   param      = "veldir2D"
   #get the Turbidity data for the corresponding time
   #time, data2D     =  Extract_matrix(df_avg_up, date_list_UTC, "veldir2D", num_range=30)
   #num_range=30
   time, r, data2D     =  Extract_matrix(df_avg, date_list, "veldir2D")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix  = PARAMSDICT[param]
   ##############################
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist,ADCP_DOWN=True, rframe= "top", **PARAMSDICT)
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
    #Plot_fig(axs[ix], time,  data, param, plot_bigger, **PARAMSDICT)
    Plot_hist(axs[ix], time, data, param, plot_bigger,  **PARAMSDICT)

#Set the axis-label for the last axis
#format='%y %m %d %H %M %S'
#axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%m-%d:%H"))
#axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%H-%M:%S"))

axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%H:%M"))

###axs[nfigs -1].xaxis.set_major_formatter(dt.DateFormatter("%H-%M"))
#myFmt       = mdates.DateFormatter('%d')

#myFmt       = mdates.DateFormatter('%d-%H')
#axs[nfigs -1].xaxis.set_major_formatter(myFmt)

#axs[nfigs -1].xaxis_date()
#fig.autofmt_xdate()

#axs[nfigs -1].xaxis.set_major_locator(HourLocator())
#axs[nfigs -1].xaxis.set_major_formatter(DateFormatter('%H:%M'))

basename_ = os.path.basename(mseedFile) 

basename = basename_.split(".mseed")[0]



#### Write on the x-axis
if(STARTDATE != ENDDATE):
    axs[nfigs -1].set_xlabel('Time from %s to %s'%(STARTDATE, ENDDATE), fontsize = 13)
    #figname = "ADCP_WAVEFORM_PLOTS_from_%s_to_%s.png"%(STARTDATE, ENDDATE)
    figname = "ADCP_%s_from_%s_to_%s.png"%(basename, STARTDATE, ENDDATE)
if(STARTDATE == ENDDATE):
    axs[nfigs -1].set_xlabel('Time of the day: %s'%(STARTDATE), fontsize = 13)
    figname   = "ADCP_%s_of_%s.png"%(basename, STARTDATE)
# Set the font size of yticks
plt.yticks(fontsize=13)
# Set the font size of xticks
plt.xticks(fontsize=14)
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
