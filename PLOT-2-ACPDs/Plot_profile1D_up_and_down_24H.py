
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
#######################################
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Affine2D
#######################
import xlrd
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import matplotlib.dates as mdates
#######################################
from matplotlib.colors import Normalize
import matplotlib.cm as cm
############################
import datetime as dtt 
import xarray as xr
import pandas as pd
from mhkit import dolfyn
from mhkit import dolfyn as dlfn
#################
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
size = 13
def InterSpline(Date_of_Day, t_in, d_in):
    #range the date for one day
    #date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods=3600)
    date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods=3600 * 24)
    #change the date_rng to numpy
    date_rng  = date_rng.to_numpy()
    #defined time step, here it's 1 second
    tstep     = 1.0
    spl       = InterpolatedUnivariateSpline(t_in, d_in)
    #Interpolation
    #t_new     = np.arange(min(t_in), max(t_in), tstep)
    #New time and data new
    t_new     = np.arange(0, date_rng.size, tstep)
    d_new     = spl(t_new)
    return (date_rng, d_new)

def hist_nbin(data):
    nbins = int(1+3.22 * np.log10(data.size))
    #nbins = np.arange(0, nbins +1, 0.15)
    nbins = np.arange(0, nbins +1, 1)
    return(nbins)


#define text size
f_size = 15
#############################
def xy_point(xmin, xmax, ymin, ymax):
    xp         = (xmin + xmax)/2
    #xp         = xmin + 0.28 * (xmax - xmin)
    #yp         = ymin + 0.8 * (ymax -ymin)
    yp         = ymin + 0.85 * (ymax -ymin)
    return(xp, yp)



def SYMETRIC_BIN(data):
    nbins = int(1 + 3.322 * np.log10(data.size)) 
    nbins_symetric = np.arange(-nbins +1, nbins +1, 0.15) 
    return(nbins)


def plot_hist(ax, data, param, **PARAMSDICT):
    ylabel, color, ix, txt = PARAMSDICT[param]
    #calcutate the nbinbs
    nbins = SYMETRIC_BIN(data)
    ax.hist(data, density= False, color = color, ec ='k', lw=3.0, alpha=1.0, rwidth = 1.0, 
            bins = nbins, histtype = 'bar',label= param.replace('_', ' '))
    pad = 6
    #get the minimum and the maximum of yaxis
    ymin, ymax    = ax.get_ylim()
    xmin, xmax    = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    xp, yp        = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ##################################
    #Add the legend
    #ax.legend(loc="upper right")
    ax.legend()
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    ax.set_xlabel('Velocity (m/s)', labelpad = pad, fontsize = 11)
#    plt.yticks(fontsize = 11)
    #number of vertical lines for grid
    locations  = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", which ='both', alpha = 0.3)



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


def read_mseed(Date_of_Day, FILE, xmlFile, bandpass):
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
    #Remove instrument response
    st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
    #st.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",water_level=60)
    #st.remove_response(inventory=inv, pre_filt=pre_filt, output="DEF",water_level=60)
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
        #Title="%s   %s  %s    %s   %s" % ( network, station, channel, stime.split(".")[0], endtime.split(".")[0])
        Title         = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
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



###################################3
def Extract_database_UTC(df,  TimeList, param):
    #Grab the time of the data frame
    time_df         = df['DateTime']
    #Transform the pandas.core.series.Series time into a list of time
    Time_ALL        = [str(t).split('.')[0] for t in time_df.to_numpy()]
    #Grab all the data corresponding to the parameter needed
    Data_ALL        = list(df[param].to_numpy())
    #Convert the time List into  to Pandas DatetimeIndex
    time_CH         = pd.to_datetime(TimeList)
    # Localize time to a specific timezone (Switzerland time)
    time_indx_utc   = time_CH.tz_localize('Europe/Zurich').tz_convert('UTC')
    #Convert DatetimeIndex to Serie
    time_series_utc = time_indx_utc.to_series()
    #Formatted the time 
    time_series_utc = time_series_utc.dt.strftime('%Y-%m-%dT%H:%M:%S')
    #Convert Time series to Timestamps
    timestamps      = pd.to_datetime(time_series_utc)
    #get the sarttime
    tsatrt          = timestamps.min()
    #Calculer le temps ecoule depuis le debut de la journee
#    tad             = (timestamps - tsatrt) % pd.Timedelta(days = 1)
    #Convert the time_series_utc in a list of string
    time_series_utc = [str(t) for t in time_series_utc]
    #id_time         = [t for t in time_data if(t in time_series_utc)]
    Select_data     = [d  for d, t in zip(Data_ALL, Time_ALL) if(t in time_series_utc)]
    return (np.asarray(time_series_utc),  np.asarray(Select_data))
    #return (np.asarray(pd.Series(tad)),  np.asarray(Select_data))


#############################################
def Extract_database(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Perform the Extraction for the correspondig Year and drop all NaN
    dYR  = [df.where(df["YR"]==float(ti.split('-')[0][-2:])).dropna(how='all') for ti in Tlist]
    #Perform the Extraction for the correspondig Month and drop all NaN
    dMTs = [dm.where(dm["MO"]==float(ti.split('-')[1])).dropna(how='all') for dm, ti  in zip(dYR,Tlist)]
    #Perform the Extraction for the correspondig Day and drop all NaN
    dDay = [dd.where(dd["DA"]==float(ti.split('-')[2])).dropna(how='all') for dd, ti  in zip(dMTs, Tlist)]
    #Now let's collect the final Data for the corresponding date for a giving paramter
    Data_ALL = np.concatenate([ds[param].to_numpy() for ds in dDay])
    #Now let's collect the final time by concatenating the times that has been collected
    Time_ALL = np.concatenate([ds['DateTime'] for ds in dDay])
    #Free Momory by deleting AYR, AMT
    del dYR, dMTs
    #id_time, sel_data = Extract_database_UTC(df,  Time_ALL, param)
    #Return the  value corresponding to the date of your choice
    return (Time_ALL, Data_ALL)
    #return (time_series_utc.to_numpy(), SEL_DATA_UTC)
    #return (id_time, sel_data)


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
        DATA2D.append(temp)
        DEPTHS.append(float(depth))

    #Trans formed the List into# 2D array
    DATA2D     = np.asarray(DATA2D)

    #Convert the TIMES into pan#das time
    TIMES      = pd.to_datetime(TIMES)

    #Check if you need everage #at the depth
    if(AVERAGE):
        data   = np.mean(DATA2D, axis = 0)
        return (TIMES, data)
    else:
        return(TIMES,DEPTHS, DATA2D)


def Extract_df_list(df,  Tlist, param, nrange=None):
#def Extract_df_list(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
    if(nrange== None):
        r        = df.range.data
    else:
        r        = df.range.data[:nrange]
    #Check if the desired velocity is the Horizontal velocity
    if(param=="velocity"):
        #print(param)
        #exit()
        #Loop over the list of the time
        try:
            P  = np.concatenate([np.nanmean(df.velds.U_mag.sel(time = ti), axis =0) for ti in Tlist])
            T  = np.concatenate([df.velds.U_mag.sel(time = ti)['time'] for ti in Tlist])
        except:
            print("The Data seems not to exist for the certain data in the list Time :  ")
            print(' , '.join(Tlist))
        #print(P)

    elif(param=="velocity_profile_up" or param=="velocity_profile_down"):
        #Loop over the list of the time
        U        = df.velds.U_mag
        Matrix2D = np.concatenate([U.sel(time = ti, range = r) for ti in Tlist], axis =1)
        try:
            #P  = np.concatenate([np.nanmean(df.velds.U_mag.sel(time = ti), axis =1) for ti in Tlist])
            ##Get the range, altitude
            P  = np.nanmean(Matrix2D ,axis =1) 
            #Get the range, altitude
            T  =  r
        except:
            print("The Data seems not to exist for the certain data in the list Time :  ")
            print(' , '.join(Tlist))
    #Check if the desired velocity is the vertical velocity
    elif(param=="vertical_vel_up"):
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df.velds.w.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df.velds.w.sel(time = ti)['time'] for ti in Tlist])

    elif(param== "vertical_vel_down"):

        P  = np.concatenate([np.nanmean(df.velds.w.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df.velds.w.sel(time = ti)['time'] for ti in Tlist])

    elif(param=="veldir1D"):
            #Loop over the list of the time
            P  = np.concatenate([np.nanmean(df.velds.U_dir.sel(time = ti), axis =0) for ti in Tlist])
            T  = np.concatenate([df.velds.U_dir.sel(time = ti)['time'] for ti in Tlist])
        
    elif(param=="avg_BS"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0) 
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])

    else:
        #Loop over the  list of the time an extract the desire parameter
        P  = np.concatenate([df[param].sel(time = ti) for ti in Tlist])
        T  = np.concatenate([df[param].sel(time = ti)['time'] for ti in Tlist])
    #############################################
    #T           =  T[: num_range]
    #P           =  P[: len(T)]
    #Return the  value corresponding to the date of your choice
    return (T, P)

#Extract 2D matrix from ADCP
#def Extract_matrix(df, Tlist,  param, nrange=25):
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
            #r        = v2d.range.data[:nrange]
            #v =====> to the NORTH
            v2d      = df.velds.v
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

        elif(param=="VERTICAL_VEL"):
            #Get the 2D Vertical velocity  
            #w =====> is the vertical component of the velocity
            #r        = w2d.range.data[:nrange]
            w2d      = df.velds.w
            w2d_sel  = np.concatenate([w2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            time     = np.concatenate([w2d.sel(time = ti)['time'] for ti in Tlist])
            #get 1D Eastern velocity
            w1d      = np.nanmean(w2d_sel.data, axis= 0)
            #Free memory
            del w2d
            #return values
            return  time, r, w1d, w2d_sel
        elif(param == "HORIZONTAL_VEL"):
            #Get the 2D Vertical velocity  
            #r        = U2d.range.data[:nrange]
            #w =====> is the vertical component of the velocity
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
            #u =====> to the EAST
            #r        = u_2d.range.data[:nrange]
            u_2d     = df.velds.u
            u2d_sel  = np.concatenate([u_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            ########################
            #Get the 2D Northen velocity  
            #v =====> to the NORTH
            v_2d     = df.velds.v
            #r        = v_2d.range.data[:nrange]
            #r        = w_2d.range.data[:nrange]
            v2d_sel  = np.concatenate([v_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            #Get the 2D Vertical velocity  
            #w =====> is the vertical component of the velocity
            w_2d     = df.velds.w
            w2d_sel  = np.concatenate([w_2d.sel(time = ti, range = r) for ti in Tlist], axis =1)
            ##########################
            time     = np.concatenate([u_2d.sel(time = ti)['time'] for ti in Tlist])
            #return values
            return   u2d_sel, v2d_sel, w2d_sel 

        elif(param=="backscatter_up"):
            alpha = 0.6
            #alpha = 0.1
            time, Matrix2D = BSC_Correct(df, Tlist, alpha, nrange, MEAN=False)
            return (time, Matrix2D)
        elif(param=="backscatter_down"):
            alpha = 2.4
            #alpha = 1.0
            time, Matrix2D = BSC_Correct(df, Tlist, alpha, nrange, MEAN=False)
            return (time, r, Matrix2D)

def Plot_Temp2D(ax, fig_object, time, data2D, depths, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
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
    y_lim_min        = min(depths)
    y_lim_max        = max(depths)
    #Set the extent for 2D plot
    extent        = [start_num , end_num,  y_lim_min, y_lim_max]
    #Make a 2D plot
    im            = ax.imshow(data2D, cmap = color, vmin = min(vmind) , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
                   interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
    
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)

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
    #ax.set_ylabel(ylabel, labelpad = None, fontsize = 11)
    ax.set_ylabel('Depth (m)', labelpad = None, fontsize = 11)
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
    ticks   = [1 * int(it) for it in depths]
    #ticks   = [it for it in depths]
    #print(ticks)
    #exit()
    #Set the inverted ticks label
    #ax.set_yticklabels(ticks)
    ax.set_yticklabels(ticks)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(base=5))
#   ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #Make the grid
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])





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
    return df_avg


def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)
    
    return vertical_locations


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



#Define function
def Plot_fig(ax, time, data, param, plot_bigger, linewidth=0.3,  **PARAMSDICT):
    #speed up figure plot
    #we set the beam angle here, 
    #beam_angle = 20.0
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    #Check the user need the plot to be bigger
    if(plot_bigger):
        ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        if(param=='avg_BS'):
            ax.plot(time, data, lw=6.0, linestyle ="-", color = color, alpha =0.5,label = param)
        else:
            ax.plot(time, data, lw=6.0, linestyle ="-", color = color, alpha =0.5,label = param.title())
    else:
        ax.plot(time, data, lw = linewidth, linestyle ="-", color = color, alpha =1.0, label = param.title())
        #ax.plot(time, data, lw=0.2, linestyle ="-", color = color, alpha =0.9, label = param.title())

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #Get the index associated with the x-axis of a specific subplot
    #x_index    = ax.get_subplotspec().get_topmost_subplotspec().colspan.start
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ############
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    ax.legend(loc="upper left")
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    plt.yticks(fontsize = 11)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
#    ax.set_xlim(time.min(), time.max())
    #Get the start and the endtime
    tstart, tend    =  start_endtime(time) 
    ax.set_xlim(tstart, tend)
    #disable axis
    ax.set_xticklabels([])
    #Extract the parameter
    file_name  = 'EXTRACTED_%s.dat'%(param)
    ##############################
    line =  ax.lines[0]
    EXTRACT_data = line.get_ydata()
    #Save the array into na file
#    np.savetxt(file_name ,PDS_MEANS , delimiter="\n", fmt="%.4f")
    np.savetxt(file_name , EXTRACT_data , delimiter="\n", fmt="%.4f")
###################################################################################

#Define function
def Plot_fig1D(ax, ranges, data, param, plot_bigger,  **PARAMSDICT):
    ##Get the parameters
    #xlabel, color, ix = PARAMSDICT[param]
    xlabel, color, ix, _ = PARAMSDICT[param]
    #ax.plot(time, data, kwargs)
    pad = 12
    #Check the user need the plot to be bigger
    if(plot_bigger):
        #ax.plot(ranges, data, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        #ax.plot(ranges, data, lw=12.0, linestyle ="-", color = color, alpha =0.5,label = param.title())
        ax.plot(data, ranges, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        ax.plot(data, ranges, lw=12.0, linestyle ="-", color = color, alpha =0.5,label = param.title())
    else:
        #ax.plot(time, data, lw=1.0, linestyle ="-", color = color, alpha =1.0, label = param.title())
        #ax.plot(ranges, data, lw=1.0, linestyle ="-", color = color, alpha =0.9, label = param.title())
        ax.plot(data, ranges,lw=1.0, linestyle ="-", color = color, alpha =0.9, label = param.title())

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #Get the index associated with the x-axis of a specific subplot
    x_index    = ax.get_subplotspec().get_topmost_subplotspec().colspan.start
    #get the coordinate of the point, where to write the text
    xp         = xmax - 0.1 * xmax
    yp         = ymax - 0.1 * ymax
    #Write the text on the figure
    ax.annotate(AphaDict[x_index] , xy=(xp, yp), xycoords='data', color='k', fontsize= 12)
        #ax.text(xp, yp, AphaDict[x_index], ha="center", verticalalignment='center', color="k", transform=ax.transAxes)

    #Invert the y-axis
    #ax.invert_yaxis()
    #ax.invert_xaxis()
    ticks         = ax.get_yticks()
    #Put ticks of xaxis on top
    ax.xaxis.tick_top()
    ############
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    ax.legend(loc="upper left")
    #put the label on top
    ax.xaxis.set_label_position('top')
    #Set labels
    ax.set_xlabel(xlabel, labelpad = pad, fontsize = 13)
    #plt.yticks(fontsize = 11)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #Set y-ticks to be negative
    depths_ticks   = ["%.1f"%(it) for it in ticks]
    ax.set_yticklabels(depths_ticks)
    #ax.set_yticklabels([-float(tick) for tick in ax.get_yticks()])
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #disable axis
#    ax.set_yticklabels([])







####################################################################################
#def Plot_1DP(ax, time, data1D, data2D_H, data2D_V, beam_angle,blank_dist, depth_adcp,  param, rframe = None, ADCP_DOWN =True,  **PARAMSDICT):
def Plot_1DP(ax, data1D, data2D_H, data2D_V, beam_angle,blank_dist, height_adcp,  param,  ADCP_DOWN =True,  **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    #Define another axis in such away that ax and ax2 share the same x-axis with different scales
#    ax2 = ax.twiny()

    Mg  = np.sqrt( data2D_H **2 + data2D_V **2)
    #Plot the paricles motion with the direction
    #ax2.quiver(data2D_H,  data2D_V, angles = 'xy', scale_units='width', scale =6, 
    param_new = "%s %s"%(param.split("_")[1], param.split("_")[2])
#    quiver = ax2.quiver(data2D_H,  data2D_V, Mg, angles = 'xy', scale_units='width', scale =6, 
#               #zorder= 2,  headwidth=3., headlength =4., cmap ='jet', color = color,   alpha =0.9, label = param)
#               zorder= 2,  headwidth=3., headlength =4., cmap ='jet', color = color,   alpha =0.9, label = param_new.upper())
#               #pivot= 'tail',   color = color, alpha =0.9, label = param)

    #Assuming rows represent depth, here we calculate the mean over columns for each row (depth)
    #depth_profile = Mg.mean(axis=1)
    depth_profile =  np.nanmean(Mg, axis = 1)
    print( Mg.shape) 
    #Step 3: Plot depth profile
    depth = np.arange(Mg.shape[0])
    print(depth_profile)
    print(len(depth_profile))
    #exit()
    #set the colorbar
    #cbar = plt.colorbar(quiver)
    #plt.colorbar(quiver, label = "Velocity (m/s)")
    #Plot the data1D this is a fake plot, just to make the time appear on the x-axis as we need
    #Use white color the curve will not appear on the particles motion plot
    #Let's now remove label on the ax2-axis so that only label on ax-axis will appear on the plot
#    ax2.axis("off")
    #ax.set_xticklabels([])

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ############
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    ax.legend(loc="upper left")
    #Set label
#    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
#    plt.yticks(fontsize = 11)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
#    ticks         = r
    #if(ADCP_OPTION):
    ax.invert_yaxis()
    if(ADCP_DOWN):
        #invert the y-axis
        #depths    = depth_adcp - np.cos(beam_angle) * ticks         
        depths     = (height_adcp+ blank_dist) - np.cos(beam_angle) * r         
        ax.invert_yaxis()
        #print(depths)
        #exit()
        #invert ticks values from positive to negative 
        #ax.set_ylabel(ylabel, labelpad = pad, loc="top", fontsize = 12)
        #move the label by -0.08 on x-axis and by 1.2 on y-axis
        #ax.yaxis.set_label_coords(-0.08 , 1.3)
        plt.yticks(fontsize = 12)
        #set the colorbar position
        #cbar.set_label("Velocity (m/s)",fontsize = 12)
        #cbar.ax.yaxis.set_label_position("right")
        #shift the label by 3.5 along x-axis and by 1.0 along the y-axis
        #cbar.ax.yaxis.set_label_coords(3.5, 1.0)
    else:
        #depths     = depth_adcp + np.cos(beam_angle) * ticks         
        depths     = (height_adcp+ blank_dist) + np.cos(beam_angle) * r         
        #ax.set_xticks([])
        #Remove the ticks marks (small vertical lines) on the x-axis
        #ax.tick_params(axis="x", length = 0, color= "white", width = 0)

    ####### Plot the  1D profile ###############
    depths_ticks   = ["%.1f"%(it) for it in depths]
    ###################### 
    depths_ticks   = sorted(["%.1f"%(it) for it in depths], reverse=True)
    #ax.plot(depth_profile, r, color='k', lw=6)
    if(ADCP_DOWN):
        #depths_ticks   = sorted(["%.1f"%(it) for it in depths], reverse=True)
        ax.plot(depth_profile, depths, color='k', lw=6, label = 'ADCP_DOWN')
    else:
        ax.plot(depth_profile, depths, color='k', lw=6, label = 'ADCP_UP')
    ####################################################
    #print(depths_ticks)
    ax.set_yticklabels(depths_ticks)
    ax.legend(loc="upper left")
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)

#    #Remove the frame on the plot
#    ax.spines["right"].set_visible(False)
#    if(rframe == "top"):
#        ax2.spines["top"].set_visible(False)
#        ax.spines["top"].set_visible(False)
#    elif(rframe == "bottom"):
#        ax2.spines["bottom"].set_visible(False)
#        ax.spines["bottom"].set_visible(False)
#    elif(rframe == "right"):
#        ax2.spines["right"].set_visible(False)
#        ax.spines["right"].set_visible(False)
#    elif(rframe =="left"):
#        ax2.spines["left"].set_visible(False)
#        ax.spines["left"].set_visible(False)
 
#    #Get the start and the endtime
#    tstart, tend    =  start_endtime(time) 
#    ax.set_xlim(tstart, tend)
    #disable axis
#    ax.set_xticklabels([])



def Plot_3DPM(ax, u2D, v2D, w2D,  param, ADCP_DOWN =True,  **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    #Define another axis in such away that ax and ax2 share the same x-axis with different scales
    #ax = plt.figure().add_subplot(projection='3d')
    ##Reshape the matrix#######################
    u3D = u2D.reshape(25,12,12)
    v3D = v2D.reshape(25,12,12)
    w3D = w2D.reshape(25,12,12)
    ##############################################
    x   = np.linspace(0, 1, 25)
   ################################
    y   = np.linspace(0, 1, 12)
   #################################
    z   = np.linspace(0, 1, 12)
    #Meshgrid the data
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    ####################################
    ax.quiver(X, Y, Z, u3D, v3D, w3D, length=0.4)

#    #Plot the paricles motion with the direction
#    quiver = ax2.quiver(data2D_H,  data2D_V, Mg, angles = 'xy', scale_units='width', scale =6, 
#               zorder= 2,  headwidth=3., headlength =4., cmap ='jet', color = color,   alpha =0.9, label = param)
#    #set the colorbar
#    plt.colorbar(quiver, label = "Velocity (m)")
    #Plot the data1D this is a fake plot, just to make the time appear on the x-axis as we need
    #Use white color the curve will not appear on the particles motion plot
    #Let's now remove label on the ax2-axis so that only label on ax-axis will appear on the plot
    ax.set_xticklabels([])

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ############
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
#    ax2.legend(loc="upper left")
    #Set label
#    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
#    plt.yticks(fontsize = 11)
#    #number of vertical lines for grid
#    locations = verical_grid()
#    ax.xaxis.set_major_locator(locations)
#    #get ticks values
#    ticks         = ax.get_yticks()
#    #if(ADCP_OPTION):
#    if(ADCP_DOWN):
#        #invert the y-axis
#        ax.invert_yaxis()
#        #invert ticks values from positive to negative 
#        ticks   = [-int(it) for it in ticks]
#    #Set the inverted ticks label
#    ax.set_yticklabels(ticks)
#    #Make the grid of the axis
#    ax.grid(visible = True, axis = "both", alpha = 0.7)
##    ax.set_xlim(time.min(), time.max())
#    #Get the start and the endtime
#    tstart, tend    =  start_endtime(time) 
#    ax.set_xlim(tstart, tend)
#    #disable axis
#    ax.set_xticklabels([])





def Plot1DCurrent_Profile(ax, ranges, data2D_H, data2D_V,  param, ADCP_DOWN =True,  **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    #Define another axis in such away that ax and ax2 share the same x-axis with different scales
    #Define 1D Horizontal velocity
    U_H1D   = np.nanmean(data2D_H, axis= 1)
    W_H1D   = np.nanmean(data2D_V, axis= 1)
    #Compute the Full velocity
    V_FULL  = np.sqrt(U_H1D **2 + W_H1D **2)
    #Select the x points that will be used for meshgrid, to determine vector postion on the plot
    xsel   = [min(V_FULL), max(V_FULL)]
    for i, vel in zip( np.arange(len(V_FULL)), sorted(V_FULL) )   :
        if(i % 5 ==0 and vel not in xsel and i != len(V_FULL) -5 ):
              xsel.append(vel)

    #Select the y points that will be used for meshgrid, to determine vector postion on the plot
    ysel   = [min(ranges), max(ranges)]
    for j, rj in zip( np.arange(len(ranges)), sorted(ranges) )   :
        if(j % 5 ==0 and rj not in ysel and j != len(ranges) -5 ):
              ysel.append(rj)
    #Performed the mesgrid for the plot
    X, Y = np.meshgrid(sorted(xsel) ,sorted(ysel))

#    ax2 = ax.twiny()

    #Plot the paricles motion with the direction
    ax.quiver(U_H1D,  W_H1D, angles = 'xy', scale_units='width', scale =1, color = color, alpha =.9, label = param)
    #Plot the data1D this is a fake plot, just to make the time appear on the x-axis as we need
    #Use white color the curve will not appear on the particles motion plot
#    ax.plot(time, data1D, color='white')
    #Let's now remove label on the ax2-axis so that only label on ax-axis will appear on the plot
#    ax2.set_xticklabels([])

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ############
#    if(float(ymax) < 0.1):
#        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText = True)) 
#        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    ax.legend(loc="upper left")
    #Set label
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    plt.yticks(fontsize = 11)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #rounded ticks
    ticks        = [round(it,2) for it in ticks]
    #if(ADCP_OPTION):
#    if(ADCP_DOWN):
#        #invert the y-axis
#        ax.invert_yaxis()
#        #invert ticks values from positive to negative 
#        #ticks   = [-1 * it for it in ticks]
#    #Set the inverted ticks label
    ax.set_yticklabels(ticks)
    #####################3
    #Apply a 90 degree rotation
#    rotation = Affine2D().rotate_deg(90)
#    #Apply the transformation to the entire axis
#    ax.transData = rotation + ax.transData
#    ax.invert_xaxis() 
#    ax.invert_yaxis()

    #Make the grid of the axis
#    ax.grid(visible = True, axis = "both", alpha = 0.7)
##    ax.set_xlim(time.min(), time.max())
#    #Get the start and the endtime
#    tstart, tend    =  start_endtime(time) 
#    ax.set_xlim(tstart, tend)
#    #disable axis
#    ax.set_xticklabels([])




########################################################################################
def Plot_fig2D(ax, fig_object, time,r, data2D, beam_angle, height_adcp,  blank_dist,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    vmind         = np.min(data2D, axis = 1) 
    vmaxd         = np.max(data2D, axis = 1)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time) 
    #set pad 
    pad = None
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
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)

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
    if('CW' in ylabel):
        #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
        if(ADCP_DOWN):
            #set the color bar on the figure 
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
            cbar_ax.yaxis.set_label_position("right")
            #shift the label by 3.5 along x-axis and by 1.0 along the y-axis
            cbar_ax.yaxis.set_label_coords(3.5, 1.1)

            #move the label by -0.08 on x-axis and by 1.3 on y-axis
            ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = 12)
            ax.yaxis.set_label_coords(-0.08 , 1.3)
        else:
            fig_object.colorbar(im, cax=cbar_ax, ticks =[0, 90, 180, 270, 350])
    #elif('CW' not in ylabel):
    else:
        if(ADCP_DOWN):
            #set the color bar on the figure 
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            cbar_ax.yaxis.set_label_position("right")
            #position of the color bar
            cbar_ax.yaxis.set_label_coords(3.5, 1.1)

            #ylabel and the postion of its positions 
            ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = 12)
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
            ax.yaxis.set_label_coords(-0.08 , 1.3)
            #plt.yticks(fontsize = 12)
        else:
            fig_object.colorbar(im, cax=cbar_ax)

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
#   ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
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
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])

##############################################

def Plot_hist(ax, time, data, param, plot_bigger,  **PARAMSDICT):
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    nbins             = hist_nbin(data)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
    ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    #xwd = 1/data.size + (1./data.size) * 0.89
    #ax.bar(time, data, width = 0.01, align ='edge')
    #compute the width as function of the data
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





def Reshape_to_3D(u, v, w):
   #set a scale to 1.   #set a scale to 1.00
   scale = 1.
   #####################
   u3D = u.reshape(25,12,12) * scale
   v3D = v.reshape(25,12,12) * scale
   w3D = w.reshape(25,12,12) * scale
   ##############################################
   x   = np.linspace(0, 1, 25)
   ################################
   y   = np.linspace(0, 1, 12)
   #################################
   z   = np.linspace(0, 1, 12)
   #Meshgrid the data
   X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
   ####################################
   return (X, Y, Z, u3D , v3D , w3D)

def Reset_axis(ax):
   ####################################################3
   ax.set_xticklabels([])
   ax.set_yticklabels([])
   ax.set_zticklabels([])
   #Set random  colors for the arrows

   ############################
   ax.set_xlabel('EAST', rotation = -15,labelpad = -10.5)
   ax.set_ylabel('NORTH', labelpad = -10.5)
   ax.set_zlabel('VERTICAL', labelpad = -10.5)
   ######set the legend ####
   ax.legend(loc = 'upper left',borderaxespad =0.005/2.0)
###############################################################3






#Open the configuration file
#with open("config2adcp.yaml") as Fym:
with open("config1D.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
PARAMSDICT     = Fig_params['PARAMSDICT']
################
PLOT2ADCP      = Fig_params['PLOT2ADCP']
#Set the SEAFLOOR DEPTH at the location of the ADCP
SEAFLOOR_DEPTH = Fig_params['SEAFLOOR_DEPTH'] # in meter
#Set the ADCP DEPTH
ADCP_DOWN_HEIGTH_FROM_LAKEBED = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']  #in meter
ADCP_UP_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']  #in meter
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFile      =   Fig_params['mseedFile']
#Seimic mseed file
Response_File  = Fig_params['Response_File']
#Grab the path of the data
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
#Grab the type of the ADCP data to plot
ADCP_DOWN      = Fig_params['ADCP_DOWN']
#Discharge file name
FILE_DISCH     = Fig_params['FILE_DISCH']
#Grab the Meteo data file
FILE_METEO     = Fig_params['FILE_METEO']
#Grab the temperature from BRB-SOLO-3
FILE_TEMPE     = Fig_params['FILE_TEMPE']
########################################
#Check the user want to make a plot bigger?
plot_bigger    = Fig_params['plot_bigger']
#Grab the starttime
STARTDATE      = Fig_params['STARTDATE']
#Plot Velocity Parameters it True or False
velocity_up    = Fig_params['velocity_up']
velocity_down  = Fig_params['velocity_down']
#Parameter to plot Horizontal velocity Direction (2D) in degrees clockwise from North: degrees_CW_from_N
veldir2D_up    = Fig_params['veldir2D_up']
veldir2D_down  = Fig_params['veldir2D_down']
#Parameter to plot Horizontal velocity Direction (1D) in degrees clockwise from North: degrees_CW_from_N
veldir1D_up    = Fig_params['veldir1D_up']
veldir1D_down  = Fig_params['veldir1D_down']
#Plot the vertical velocity Parameters it True or False
vertical_vel_up = Fig_params['vertical_vel_up']
vertical_vel_down = Fig_params['vertical_vel_down']

################Plot 1D Velocity Profile #######################
velocity_profile_down = Fig_params['velocity_profile_down']
velocity_profile_up   = Fig_params['velocity_profile_up']
#Plot Pressure Parameters it True or False
pressure_up     = Fig_params['pressure_up']
pressure_down   = Fig_params['pressure_down']
###################################
#Plot the depth change of ADCP
depth_adcp_up   = Fig_params['depth_adcp_up']
depth_adcp_down = Fig_params['depth_adcp_down']

#Plot particle motions in the water column
pm_adcp_up     = Fig_params['pm_adcp_up'] 
pm_adcp_down   = Fig_params['pm_adcp_down'] 

#############################################
hist_adcp_down = Fig_params['hist_adcp_down'] 
hist_adcp_up   = Fig_params['hist_adcp_up'] 

############################################
pm3D_adcp_up   = Fig_params['pm3D_adcp_up']
pm3D_adcp_down = Fig_params['pm3D_adcp_down']
##################################################
###########Surface tracking #########################
surface_tracking_up   = Fig_params['surface_tracking_up']
surface_tracking_down = Fig_params['surface_tracking_down']

#1D Current profile 
Current_Profile_up   = Fig_params['Current_Profile_up']
Current_Profile_down = Fig_params['Current_Profile_down']

##Grab the Velocity component to visualize
temperature_up    = Fig_params['temperature_up']
temperature_down    = Fig_params['temperature_down']
############################################
Temperature2D  = Fig_params['Temperature2D']
##Grab the wind speed paramter to plot
wind_speed     = Fig_params['wind_speed']
##Grab the wind speed direction  paramter to plot
wind_direction     = Fig_params['wind_direction']
#Get the precipation parameter
precipitation  = Fig_params['precipitation']
#Grab the discharge parameter
discharge      = Fig_params['discharge']
#Plot the Average Backscatter?
avg_BS_up        = Fig_params['avg_BS_up']
avg_BS_down        = Fig_params['avg_BS_down']
##Plot the backscatter of Particules ###
backscatter_up = Fig_params['backscatter_up']
backscatter_down = Fig_params['backscatter_down']
#Plot the 2D velocity
velocity2D_up   = Fig_params['velocity2D_up']
velocity2D_down = Fig_params['velocity2D_down']
# Plot the seismic Waveform 
waveform       = Fig_params['waveform']
#Plot the seismic Envelope
envelope       = Fig_params['envelope']
#Plot seismic spectrogram
spectrogram    = Fig_params['spectrogram']
#Grab the bandpass
bandpass       = Fig_params['bandpass']
bandpass_spec  = Fig_params['bandpass_spec']
#get the nbin
n_bin          = Fig_params['n_bin']
#Get the number of range 2D plot 
num_range      = Fig_params['num_range']

#Grab the time list
try:
    date_all_UTC         = pd.date_range(start= STARTDATE, end = STARTDATE, tz = 'UTC')
    date_all_CH          = pd.date_range(start= STARTDATE, end = STARTDATE, tz = 'Europe/Zurich')
    #Convert the Switzerland Time into UTC time
    #date_all_CH2UTC      = date_all_CH.tz_convert('UTC')
except:
     print("Checked the date entring %s"%(STARTDATE))
     exit()
#Create a time list for the started and the
date_list_UTC = [str(it).split()[0] for it in date_all_UTC]
date_list_CH  = [str(it).split()[0] for it in date_all_CH]
#Convert the Switzerland Time into UTC time
################################# Provide the files ###################


#Set the number of subfigures to Zero
nfigs = 0
#Loop over the dictionary parameters
for i, key_param in zip(range(len(PARAMSDICT)), PARAMSDICT):
    #Check if the parameter is plotable
    #print(key_param)
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
        #Add the corresponding anotation letter
        PARAMSDICT[key].append(AphaDict[i])

#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
nfigs  = len(PARAMSDICT)


#Read the ADCP file
try:
    if(PLOT2ADCP):
        #read ADCP looking downward
        df_down      = dlfn.read(ADCP_FILE_NAME_DOWN) 
        #read ADCP looking upward
        df_up        = dlfn.read(ADCP_FILE_NAME_UP)
        EXT          = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1] 
        #########################################
    #elif(ADCP_DOWN):
    #    df      = dlfn.read(ADCP_FILE_NAME_DOWN) 
    #    EXT     = os.path.basename(ADCP_FILE_NAME_DOWN).split("_")[1] 
    #else:
    #    df      = dlfn.read(ADCP_FILE_NAME_UP)
    #    EXT     = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1] 
############################################
except:
    print("Provide the two  ADCP-Files that need to be plot") 
    exit()


######Perform the averaging#
df_avg_down     = Average(df_down, n_bin)
df_avg_up       = Average(df_up, n_bin)
#get the beam angle
beam_angle_up   = df_up.beam_angle
beam_angle_down = df_down.beam_angle
#get the blink distance
blank_dist_up   = df_up.blank_dist
blank_dist_down = df_down.blank_dist
############################################
range_up        = df_up.range.data
range_down      = df_down.range.data
#Create the figure


#Read the Discharge  file
d_disch= pd.read_fwf(FILE_DISCH, delimiter='\t')
#Add the of DateTime
d_disch= Add_Column(d_disch)
#print(d_disch) 
#read the Meteo database
d_mteo = pd.read_csv(FILE_METEO, delimiter ='\t')
#Add the of DateTime
d_mteo = Add_Column(d_mteo)
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


if(pm3D_adcp_up == False or pm3D_adcp_down == False):
    #Create  a normal 2D figure
    fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)


if(velocity_profile_up):
    #Create the figure
     #Grab the figure size
     fig, axs  = plt.subplots(1, nfigs, sharex = False, figsize = fig_size)
#

##Create the figure
#if(pm3D_adcp_up):
#    # set up a figure twice as wide as it is tall
#    fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#
#    axs[0] = fig.add_subplot(1, 3, 3, projection='3d')
#    axs[1] = fig.add_subplot(2, 3, 3, projection='3d')
#    ##########################################3
#    # set up a figure twice as wide as it is tall
#else:
#    fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
###############################################################################
    #fig, axs  = plt.subplots(1, nfigs, sharex = False, figsize = fig_size)




#get the basename
basename_ = os.path.basename(mseedFile) 
basename  = basename_.split(".mseed")[0]
#######################################


####################################################
if(envelope):
    #set the parameter
    param      = "envelope"
    #if(spectrogram):
    lw = 1.2
    #get the seismogram from mseedFile
    #time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass, Pressure = True)
    time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Set the Title
    basename_new = basename.split('_')
    label        = '%s   %s '%(basename_new[0], basename_new[1])
    #label        = '%s   %s             BP: %s s'%(basename_new[0], basename_new[1], bandpass)
    axs[ix].set_title(label, loc='center', pad=None)
    #Plot the figure by calling the plotting function, PlotEnvelop
    #Plot_fig(axs[ix], time, tr.data, param, False, lw, **PARAMSDICT)
    #PlotEnvelop(axs[ix], STARTDATE, param, tr, linewidth=1.2,  **PARAMSDICT):
    PlotEnvelop(axs[ix], STARTDATE, param, tr, lw,  **PARAMSDICT)

#################################################################
if(waveform):
    #set the parameter
    param      = "waveform"
    if(spectrogram):
        bandpass   = bandpass_spec
        #bandpass   = bandpass
        lw = 0.3
    else:
        lw = 0.8
    #get the seismogram from mseedFile
    time, tr, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Set the Title
    basename_new = basename.split('_')
    #label        = '%s   %s             BP: %s s'%(basename_new[0], basename_new[1], bandpass)
    label        = '%s       %s    '%(basename_new[0], basename_new[1])
    if(envelope==False):
        axs[ix].set_title(label, loc='center', pad=None)
    #Plot the figure by calling the plotting function, plot_twinx
    data  =  tr.data * 1e+6
    #time = time.astype("datetime64[ns]")
    Plot_fig(axs[ix], time, data, param, False, lw, **PARAMSDICT)
    #Plot_fig(axs[ix], time_new, data, param, False, **PARAMSDICT)


if(spectrogram):
    #set the parameter
    param            = "spectrogram"
    #get the seismogram from mseedFile
    time, tr, Title  = read_mseed(STARTDATE, mseedFile, Response_File, bandpass_spec)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    plot_spectrogram(axs[ix], fig, tr, STARTDATE, **PARAMSDICT)


#Plot the Discharge
if(discharge):
    #set the parameter
    param         = "discharge"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_disch,  date_list_CH, 'Discharge')
    ####################3
    #np.savetxt('Discharge.dat' , data , delimiter="\n", fmt="%.4f")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

#Doubleched if we need to plot the Velocity
if(velocity_up):
   #set the parameter
   param      = "velocity_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "velocity")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function
if(velocity_down):
   #set the parameter
   param      = "velocity_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, "velocity")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function



if(hist_adcp_up):
   #set the parameter
   param      = "hist_adcp_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "velocity")
   #time, data = Extract_df_list(df_avg_up , date_list_UTC, 'velocity', nrange=25)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #plot_hist(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   plot_hist(axs[ix], data, param, **PARAMSDICT)

if(hist_adcp_down):
   #set the parameter
   param      = "hist_adcp_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, "velocity")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   plot_hist(axs[ix], data, param, **PARAMSDICT)
   #plot_hist(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(vertical_vel_up):
       #set the parameter
       param      = "vertical_vel_up"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg_up , date_list_UTC, param)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix, _ = PARAMSDICT[param]
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(vertical_vel_down):
       #set the parameter
       param      = "vertical_vel_down"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg_down , date_list_UTC, param)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix, _ = PARAMSDICT[param]
       data = 1 * data
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
####################################################################
if(veldir1D_up):
   #set the parameter
   param      = "veldir1D_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function

if(veldir1D_down):
   #set the parameter
   param      = "veldir1D_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, param)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function


#Plot particle motions in the water column
if(pm_adcp_up):
   param      = "pm_adcp_up"
   #get the Turbidity data for the corresponding time
   #time,r, u1d,  u2d     =  Extract_matrix(df_avg_up, date_list_UTC, "EAST_VEL")
   time, r, u1d,  u2d     =  Extract_matrix(df_avg_up, date_list_UTC, "HORIZONTAL_VEL")
   time2,r, w1d, w2d     =  Extract_matrix(df_avg_up, date_list_UTC, "VERTICAL_VEL")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
    #Call the function to plot particule motion
    #remove a bottom frame
   rframe = "bottom"
   Plot_PM(axs[ix], time, u1d, u2d, w2d, beam_angle_up, blank_dist_up,  ADCP_UP_HEIGTH_FROM_LAKEBED,  param, rframe, ADCP_DOWN =False,  **PARAMSDICT)
   #Plot_PM(axs[ix], time, u1d, u2d, w2d, beam_angle_up, blank_dist_up,  ADCP_UP_HEIGTH_FROM_LAKEBED,  param, ADCP_DOWN =False,  **PARAMSDICT)

if(pm_adcp_down):
   param      = "pm_adcp_down"
   #get the Turbidity data for the corresponding time
   #time,r, u1d,  u2d     =  Extract_matrix(df_avg_down, date_list_UTC, "EAST_VEL")
   time, r, u1d,  u2d     =  Extract_matrix(df_avg_down, date_list_UTC, "HORIZONTAL_VEL")
   time2, r,  w1d, w2d     =  Extract_matrix(df_avg_down, date_list_UTC, "VERTICAL_VEL")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   # change the sign of the vertical velocity
   w2d = -1 * w2d
   #Call the function to plot particule motion
   #remove frame
   rframe = "top"
   Plot_PM(axs[ix], time, u1d, u2d, w2d, beam_angle_down, blank_dist_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  param, rframe, ADCP_DOWN =True,  **PARAMSDICT)
   #Plot_PM(axs[ix], time, u1d, u2d, w2d, beam_angle_down, blank_dist_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  param, ADCP_DOWN =True,  **PARAMSDICT)
######################################################################
if(Current_Profile_up):
   param      = "Current_Profile_up"
   #get the Turbidity data for the corresponding time
   #time, u1d,  u2d     =  Extract_matrix(df_avg_up, date_list_UTC, "EAST_VEL")
   time, r, U1d,  U2d     =  Extract_matrix(df_avg_up, date_list_UTC, "HORIZONTAL_VEL")
   time2, r, w1d, w2d     =  Extract_matrix(df_avg_up, date_list_UTC, "VERTICAL_VEL")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
    #Call the function to plot particule motion
   Plot1DCurrent_Profile(axs[ix], r, U2d, w2d,  param, ADCP_DOWN =False,  **PARAMSDICT)
##########################################################################
if(Current_Profile_down):
   param      = "Current_Profile_down"
   #get the Turbidity data for the corresponding time
   time,r, U1d,  U2d     =  Extract_matrix(df_avg_down, date_list_UTC, "HORIZONTAL_VEL")
   time2,r, w1d, w2d     =  Extract_matrix(df_avg_down, date_list_UTC, "VERTICAL_VEL")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
    #Call the function to plot particule motion
   Plot1DCurrent_Profile(axs[ix], r, U2d, w2d,  param, ADCP_DOWN =True,  **PARAMSDICT)

####################################################

if(velocity_profile_up):
   ##################
   param = "velocity_profile_up" 
   #get the vertical velocity profile data for the corresponding time
   #ranges, data = Extract_df_list(df_avg_up , date_list_UTC, param)
   time, r, u1d,  u2d     =  Extract_matrix(df_avg_down, date_list_UTC, "HORIZONTAL_VEL")
   time2, r,  w1d, w2d     =  Extract_matrix(df_avg_down, date_list_UTC, "VERTICAL_VEL")
   #print(data)
   #exit()
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   # change the sign of the vertical velocity
   #w2d = -1 * w2d
   #Call the function to plot particule motion
   #Compute the height from the lakebed
   H_up  = (ADCP_UP_HEIGTH_FROM_LAKEBED + blank_dist_up) + np.cos(beam_angle_up) * r
   #ADCP_DOWN_HEIGTH_FROM_LAKEBED
   #Plot the figure by calling the plotting function
   #Plot_fig1D(axs[ix], ranges, data, param, plot_bigger,  **PARAMSDICT)
   Plot_1DP(axs[ix], u1d, u2d, w2d, beam_angle_up, blank_dist_up, ADCP_UP_HEIGTH_FROM_LAKEBED,  param, ADCP_DOWN =False,  **PARAMSDICT)
   #Plot_fig1D(axs[ix], H_up, data, param, plot_bigger,  **PARAMSDICT)

if(velocity_profile_down):
   ##################
   param = "velocity_profile_down" 
   #get the vertical velocity profile data for the corresponding time
   #ranges, data = Extract_df_list(df_avg_down , date_list_UTC, param)
   time, r, u1d,  u2d     =  Extract_matrix(df_avg_down, date_list_UTC, "HORIZONTAL_VEL")
   time2, r,  w1d, w2d     =  Extract_matrix(df_avg_down, date_list_UTC, "VERTICAL_VEL")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #################
   #data    = 100 * data
   # change the sign of the vertical velocity
   w2d = -1 * w2d
   #Compute the height from the lakebed
   H_down  = (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist_down) + np.cos(beam_angle_down) * r
   #print(ranges)
   #print(H_down)
   #exit()
   #Plot the figure by calling the plotting function
   #Plot_fig1D(axs[ix], H_down, data, param, plot_bigger,  **PARAMSDICT)
   Plot_1DP(axs[ix], u1d, u2d, w2d, beam_angle_down, blank_dist_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED,  param, ADCP_DOWN =True,  **PARAMSDICT)
   #Plot_1DP(axs[ix], u1d, u2d, w2d, beam_angle_up, blank_dist_up, depth_adcp_up,  param, ADCP_DOWN =False,  **PARAMSDICT)
   #Plot_fig1D(axs[ix], ranges, data, param, plot_bigger,  **PARAMSDICT)



if(veldir2D_up):
   #set the parameter
   param      = "veldir2D_up"
   #get the Turbidity data for the corresponding time
   #time, data2D     =  Extract_matrix(df_avg_up, date_list_UTC, "veldir2D", num_range=30)
   #num_range=30
   time, r, data2D     =  Extract_matrix(df_avg_up, date_list_UTC, "veldir2D")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   ##############################
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_fig2D(axs[ix],fig, time, r , data2D, beam_angle_up,ADCP_UP_HEIGTH_FROM_LAKEBED,blank_dist_up,ADCP_DOWN = False, rframe= "bottom", **PARAMSDICT)
######################
if(veldir2D_down):
   #set the parameter
   param      = "veldir2D_down"
   #get the Turbidity data for the corresponding time
   #time, r, data2D     =  Extract_matrix(df_avg_down, date_list_UTC,  "veldir2D", num_range)
   time, r, data2D     =  Extract_matrix(df_avg_down, date_list_UTC,  "veldir2D")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function of 2D Matrix
   Plot_fig2D(axs[ix],fig, time,r, data2D, beam_angle_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist_down,ADCP_DOWN = True, rframe= "top", **PARAMSDICT)
#######################################################################

   #set the parameter
if(depth_adcp_up):
   param      = "depth_adcp_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "depth")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   data = SEAFLOOR_DEPTH - data
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(depth_adcp_down):
   #set the parameter
   param      = "depth_adcp_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, "depth")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   data = SEAFLOOR_DEPTH - data
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)


if(surface_tracking_up):
   param      = "surface_tracking_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "depth")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #################################################################
if(surface_tracking_down):
   param      = "surface_tracking_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, "depth")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
   #Plot the figure by calling the plotting function

##Doubleched if we need to plot the Temperature
#if(temperature):
#   #set the parameter
#   param            = "temperature"
#   #De-comment this if you want to plot ADCP temparature
#   #get the horizontal velocity data for the corresponding time
##   time, data   = Extract_df_list(df_avg , date_list_UTC, "temp")
#   time, data       = Extract_RBR_SOLO_TEMP(data_temp,  date_list_UTC, AVERAGE=True)
#   #get the parameter, the index of the corresponding axis, and the color
#   _ , color, ix, _ = PARAMSDICT[param]
#   #Plot the figure by calling the plotting function
#   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)


if(temperature_up):
   #set the parameter
   param            = "temperature_up"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_up , date_list_UTC, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(temperature_up):
   #set the parameter
   param            = "temperature_down"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_down , date_list_UTC, "temp")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(Temperature2D):
   param            = "Temperature2D"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
#   time, data   = Extract_df_list(df_avg , date_list_UTC, "temp")
   time, depths, data2D       = Extract_RBR_SOLO_TEMP(data_temp,  date_list_UTC, AVERAGE=False)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
    #Plot the figure by calling the plotting function of 2D Matrix
   #L = fig.axes
   #print(fig.axes) 
   #print(L[0]) 
   #exit()
   #Plot_Temp2D(axs[ix],fig, time, data2D, depths, **PARAMSDICT)
   Plot_Temp2D(axs,fig, time, data2D, depths, **PARAMSDICT)
##Plot the Pressure
if(pressure_up):
   #set the pad_echr
   param         = "pressure_up"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg_up , date_list_UTC, "pressure")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(pressure_down):
   #set the pad_echr
   param         = "pressure_down"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg_down , date_list_UTC, "pressure")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
##Plot the backscatter
if(backscatter_up):
    #set the parameter
    param      = "backscatter_up"
    #get the Turbidity data for the corresponding time
    time,r, data2D     =  Extract_matrix(df_avg_up, date_list_UTC, param, num_range)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix],fig, time, data2D, num_range,  ADCP_DOWN = False, **PARAMSDICT)

if(backscatter_down):
    #set the parameter
    param      = "backscatter_down"
    #get the Turbidity data for the corresponding time
    time,r, data2D     =  Extract_matrix(df_avg_down, date_list_UTC, param, num_range)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    Plot_fig2D(axs[ix],fig, time, data2D, num_range,  ADCP_DOWN = True, **PARAMSDICT)
#Plot the 2D velocity
if(velocity2D_up):
    #set the parameter
    param      = "velocity2D_up"
    #get the Turbidity data for the corresponding time
    #time, r, data2D     =  Extract_matrix(df_avg_up, date_list_UTC,  "velocity2D", num_range)
    time, r, data2D     =  Extract_matrix(df_avg_up, date_list_UTC,  "velocity2D")
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig, time, data2D, num_range, ADCP_DOWN = False, **PARAMSDICT)
    Plot_fig2D(axs[ix],fig, time,r, data2D, beam_angle_up, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist_up,ADCP_DOWN = False, rframe= "bottom", **PARAMSDICT)
#Plot the 2D velocity
if(velocity2D_down):
    #set the parameter
    param      = "velocity2D_down"
    #get the Turbidity data for the corresponding time
    #time,r, data2D     =  Extract_matrix(df_avg_down, date_list_UTC,  "velocity2D", num_range)
    time,r, data2D     =  Extract_matrix(df_avg_down, date_list_UTC,  "velocity2D")
    #print("****" *40)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function of 2D Matrix
    #Plot_fig2D(axs[ix],fig, time, data2D, num_range, ADCP_DOWN = True, **PARAMSDICT)
    Plot_fig2D(axs[ix],fig, time,r, data2D, beam_angle_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, blank_dist_down,ADCP_DOWN = True, rframe= "top", **PARAMSDICT)
##Plot the Turbidity
if(avg_BS_up):
    #set the parameter
    param      = "avg_BS_up"
    #get the Turbidity data for the corresponding time
    time, data = Extract_df_list(df_avg_up,  date_list_UTC, "avg_BS")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
##Plot the Turbidity
if(avg_BS_down):
    #set the parameter
    param      = "avg_BS_down"
    #get the Turbidity data for the corresponding time
    time, data = Extract_df_list(df_avg_down,  date_list_UTC, "avg_BS")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(wind_speed):
    #set the parameter
    param         = "wind_speed"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_CH, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)

if(wind_direction):
    #set the parameter
    param         = "wind_direction"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_CH, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], time, data, param, plot_bigger, **PARAMSDICT)
#Plot the precipitation
if(precipitation):
    #set the parameter
    param         = "precipitation"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_CH, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _= PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    # convert hPa to bar;  1hPa =   0.001 bar
    #Plot_fig(axs[ix], time,  data, param, plot_bigger, **PARAMSDICT)
    Plot_hist(axs[ix], time, data, param, plot_bigger,  **PARAMSDICT)

################################################################################
#if(pm3D_adcp_up == True and pm3D_adcp_down = True):
if(pm3D_adcp_up):
   param    = "pm3D_adcp_down"
   #fig      = plt.figure(figsize = fig_size)
   fig3D      = plt.figure(figsize = (8,12))
   ax2      = fig3D.add_subplot(121, projection = '3d')
   ax1      = fig3D.add_subplot(221, projection = '3d')
   ##########################################
   u,  v, w =  Extract_matrix(df_avg_up, date_list_UTC, "3Comps")
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   ############Reshape the data into 3D################
   X, Y, Z, u3P, v3P, w3P = Reshape_to_3D(u, v, w)
   ##########
   #ax1.quiver(X, Y, Z, u3D *scale, v3D *scale, w3D * scale, length = 0.8, alpha = .8, lw=0.5, label = 'Upward ADCP')
   q1= ax1.quiver(X, Y, Z, u3P, v3P , w3P, color ='b',   alpha = .8, lw=0.5, label = 'Upward ADCP')
   #q1= ax1.quiver(X, Y, Z, u3D *scale, v3D *scale, w3D * scale, cmap = 'jet')
   #### Reset Axis ###################
   Reset_axis(ax1)
   ####################################################3
   if(pm3D_adcp_down):
        param    = "pm3D_adcp_down"
        ###########################################
        ud,  vd, wd =  Extract_matrix(df_avg_down, date_list_UTC, "3Comps")
        #get the parameter, the index of the corresponding axis, and the color
        _ , color, ix, _ = PARAMSDICT[param]

        ############Reshape the data into 3D################
        XX, YY, ZZ, u3d, v3d, w3d = Reshape_to_3D(ud, vd, wd)
        ####################################
        #ax2.quiver(X, Y, Z, u3d * scale, v3d * scale, w3d *scale, length=0.8 ,cmap = 'Reds', lw=0.8, alpha = .8, label='Downward ADCP')
        ax2.quiver(XX, YY, ZZ, u3d , v3d, w3d, color ='b', lw=0.8, alpha = .8, label='Downward ADCP')
        #fig.subplots_adjust(hspace = .8)
        #ax1.legend(loc = 'upper left',borderaxespad =0.005 )
        #ax2.legend(loc = 'upper right' ,borderaxespad =0.005 )
        Reset_axis(ax2)
        ####################################################
#   ax2.text2D(0.05, 0.95, 'Upward Looking ADCP', transform = ax2.transAxes, ha = 'center', va = 'center')
#   ax1.text2D(0.05, -1.0, 'Downward Looking ADCP', transform = ax2.transAxes, ha = 'center', va = 'center')
   fig3D.subplots_adjust(top = 0.9, hspace = .5)
   plt.savefig('Sample.png', bbox_inches = 'tight', dpi = 300)

##Set the ticks on the last x-axis
#if(pm3D_adcp_up == False):
if(Temperature2D):
    axs.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs.set_xlabel('Depth (m) ', fontsize = 13)
    figname   = "ADCP-%s_%s_of_%s.png"%(EXT, basename, STARTDATE)
    fig.savefig(figname, bbox_inches = 'tight', dpi = 300)

elif(velocity_profile_up == True or velocity_profile_down == True or hist_adcp_down ==True):
    figname   = "ADCP-%s_%s_of_%s.png"%(EXT, basename, STARTDATE)
    #Set the font size of yticks
    plt.yticks(fontsize=13)
    # Set the font size of xticks
    plt.xticks(fontsize=14)
    if(hist_adcp_down ==True):
        #plt.subplots_adjust(hspace = fig_space)
        plt.subplots_adjust(hspace = fig_space + 3 * fig_space)
    ##Space between the subplots
    elif(nfigs <= 5):
        #plt.subplots_adjust(hspace = fig_space)
        plt.subplots_adjust(hspace = fig_space)
    elif(nfigs > 5):
        fig_space = fig_space + (1.0/100) * nfigs
        plt.subplots_adjust(hspace = fig_space)
    
    #plt.subplots_adjust(hspace = fig_space)
    #Aligned all the ylabels of the figure
    fig.align_ylabels()
    #Save the figure
    fig.savefig(figname, bbox_inches = 'tight', dpi = 300)

elif(pm3D_adcp_up == False or pm3D_adcp_down == False):
    if(Current_Profile_up== False and Current_Profile_down==False):
        if(nfigs >= 1):
            axs[nfigs -1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            #axs[nfigs -1].xaxis.set_major_formatter(mdates.DateFormatter("%d/%H"))
            #Text to Write on the x-axis
            axs[nfigs -1].set_xlabel('Time (hour: minute) on %s'%(STARTDATE), fontsize = 13)
        else:
            #axs[nfigs].xaxis.set_major_formatter(mdates.DateFormatter("%d/%H"))
            axs[nfigs].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            #Text to Write on the x-axis
            axs[nfigs].set_xlabel('Time (hour: minute) on %s'%(STARTDATE), fontsize = 13)
    else:
        axs[nfigs -1].set_xlabel('range (m) ', fontsize = 13)
    #Set the figure name
    if(spectrogram):
        #figname   = "ADCP-%s_%s_of_%s.pdf"%('SEPC', basename, STARTDATE)
        figname   = "ADCP-%s_%s_of_%s.png"%('SEPC', basename, STARTDATE)
    else:
        #figname   = "ADCP-%s_%s_of_%s.pdf"%(EXT, basename, STARTDATE)
        figname   = "ADCP-%s_%s_of_%s.png"%(EXT, basename, STARTDATE)
    #Set the font size of yticks
    plt.yticks(fontsize=13)
    # Set the font size of xticks
    plt.xticks(fontsize=14)
    #Make a rotation of xticks
    #plt.xticks(rotation=45)
    
    ##Space between the subplots
    if(nfigs <= 5):
        plt.subplots_adjust(hspace = fig_space)
    elif(nfigs > 5):
        fig_space = fig_space + (1.0/100) * nfigs
        plt.subplots_adjust(hspace = fig_space)
    
    #plt.subplots_adjust(hspace = fig_space)
    #Aligned all the ylabels of the figure
    fig.align_ylabels()
    #Save the figure
    fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
    #fig.savefig(figname, bbox_inches = 'tight', dpi = 510)
    
#    fig.savefig(figname)
