
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
from scipy.signal import detrend
#################################################
from matplotlib.ticker import FormatStrFormatter, NullFormatter
#from matplotlib.ticker import NullFormatter
##############################################
import matplotlib.dates as mdates
#######################################
from matplotlib.colors import Normalize
import matplotlib.cm as cm
############################
import datetime as dtt 
import xarray as xr
import pandas as pd
#from mhkit import dolfyn
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



#Read the mseed file
def read_mseed(Date_of_Day, mseedFile, xmlFile, bandpass,resampling=True, Pressure = None):
    #get the stream by reading the mseed file
    st      = read(mseedFile)
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
        #################################
        #sampling_rate = tr.stats.sampling_rate
        if(resampling):
            #The sampling rate is 1 second
            #Resampling the data
            tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
            #get the data
            data = tr.data
            date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= tr.data.size)
            #change the date_rng to numpy
            time      = date_rng.to_numpy()
        else:
            data          = tr.data
            time          = tr.times()
        #return(time, tr, Title)
    return(time, data, Title)




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
    num_vertical = 3
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
def Plot_fig(ax, freqs, data, param, plot_loglog, linewidth=1.3,  **PARAMSDICT):
    #speed up figure plot
    #we set the beam angle here, 
    #beam_angle = 20.0
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    if("velocity_up" in param):
        param = "UPWARD CURRENT"
        #param = "CURRENT"
    elif("velocity_down" in param):
        #param = "DOWNWARD CURRENT"
        param = "CURRENT"
    elif("wind" in param):
        param = "WIND"
    #Check the user need the plot to be bigger
    if(plot_loglog):
        #ax.plot(freqs, data, lw=1.0, linestyle ="-", color = color, alpha =1.0)
        ax.loglog(freqs, data, lw=linewidth,  color = color, alpha =1.0)
        if(param=='avg_BS'):
            ax.loglog(freqs, data, lw=linewidth,  color = color, alpha =0.5,label = param)
        else:
            ax.plot(freqs, data, lw=linewidth, color = color, alpha =0.5,label = param.title())
            ax.loglog(freqs, data, lw=6.0, linestyle ="-", color = color, alpha =0.5,label = param)
    else:
        ax.plot(freqs, data, lw = linewidth, color = color, alpha =1.0, label = param)
        #ax.loglog(freqs, data/max(data), lw = linewidth, color = color, alpha =1.0, label = param)
        #ax.set_xscale('log')

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
    ax.legend(loc="upper right")
    #Set label
    #ax.set_ylabel(ylabel, labelpad = pad, fontsize = 11)
    #ax.set_ylabel("Power (dB/Hz)**2", labelpad = pad, fontsize = 12)
    #ax.set_ylabel(r'Power $(\mathrm{dB/Hz})^2$', labelpad = pad, fontsize = 14)
    ax.set_ylabel(r'PSD $(\mathrm{m^2s^{-2}/Hz})$', labelpad = pad, fontsize = 14)
    #fontsize of the yticks
    ax.tick_params(axis='both', labelsize=14)
    #number of vertical lines for grid
    locations = verical_grid()
#    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #plot the log scale
    #ax.set_ylim(1e-5,1.25)
    #if("wind" in param or "WIND" in param):
    #    ax.set_ylim(1e-1 * 8, 1e+5)
    #else:
    #    ax.set_ylim(1e-4, 1e+1 * 2)

#    ax.set_xlim(0.0000, 0.05)
    #ax.set_xlim(0.05, 0.5)
    #ax.set_xlim(0.03, 0.5)
    #ax.set_xlim(0.045, 0.5)
    #ax.set_xlim(0.045, 0.5)
    #ax.invert_yaxis()
    if(plot_loglog):
        #linear scale
        if("wind" in param or "WIND" in param):
            #ax.set_ylim(1e+1 * 0.08,  1e+3)
            #ax.set_ylim(1e+0 * 0.6,  1e+3)
            #ax.set_ylim(1e+0 * 0.6,  1e+4 )
            #ax.set_ylim(1e+0 * 0.6,  1e+3 * 2)
            #ax.set_ylim(1e+0 * 0.6,  1e+4 * 2)
            ax.set_ylim(1e+0 * 0.6,  1e+5 * 1.5)
        else:
            #ax.set_ylim(1e-4 *0.2, 1e+1)
            #ax.set_ylim(1e-4 *0.8, 3)
            ax.set_ylim(1e-4 *0.8, 15)
            #ax.set_ylim(1e-4 *0.8, 8)
            #ax.set_ylim(1e-4 *0.8, 2)
            #ax.set_ylim(1e-7, 0.1 *2)
#        #disable axis
        if(ix < len(PARAMSDICT) -1):
            #ax.xaxis.set_visible(False)
            #ax.xaxis.set_ticklabels([])
            # Remove x-axis labels using NullFormatter
            ax.xaxis.set_major_formatter(NullFormatter())
            ##Set ticks and labels 
        #else:
        #    #Set custom x-axis ticks and labels
        #    custom_ticks = [0.1, 0.2, 0.3, 0.4]  # Positions in terms of the base-10 log scale
        #    custom_labels = [r"$10^{-1}$", r"$2 \times 10^{-1}$", r"$3 \times 10^{-1}$", r"$4 \times 10^{-1}$"]
        #    ax.set_xticks(custom_ticks)
        #    ax.set_xticklabels(custom_labels)
        #############################################################



def Plot_waveforms_freq(ax, freqs, data, param, plot_loglog=False, linewidth=0.3,  **PARAMSDICT):
    #speed up figure plot
    print(max(data))
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix, txt = PARAMSDICT[param]
    pad = None
    #Check the user need the plot to be bigger
    if(plot_loglog):
        ax.loglog(freqs, data, lw = linewidth, linestyle ="-", color = color, alpha =1.0, label = param)
    else:
        #ax.loglog(freqs, data, lw = linewidth, linestyle ="-", color = color, alpha =1.0, label = param)
        ax.plot(freqs, data, lw = linewidth, linestyle ="-", color = color, alpha =1.0, label = param)
        ax.set_xscale('log')

    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim() 
    xmin, xmax = ax.get_xlim() 
    #Get the index associated with the x-axis of a specific subplot
    #x_index    = ax.get_subplotspec().get_topmost_subplotspec().colspan.start
    #get the coordinate of the point, where to write the text
    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    ##Write the text on the figure
    #ax.annotate(txt , xy=(xp, yp), fontsize = f_size)
    ############
    if(float(ymax) < 0.01):
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Add the legend
    #ax.legend(loc="upper right")
    #Set label
    #ax.set_ylabel(r'PSD $(\mathrm{m^2s^{-2}/Hz})$', labelpad = pad, fontsize = 14)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = 14)
    ax.set_xlabel('Frequency (Hz) ', fontsize = 14)
    #fontsize of the yticks
    ax.tick_params(axis='both', labelsize=14)
    #number of vertical lines for grid
    locations = verical_grid()
#    ax.xaxis.set_major_locator(locations)
    #Make the grid of the axis
#    ax.grid(visible = True, axis = "both", alpha = 0.7)
    #plot the log scale
    #ax.set_xlim(0.001, 1.0)    #good for loglog scale
    #linear scale
    if("waveform"==param):
        if(plot_loglog):
            ax.set_ylim(1e-11,  1e-1)  #for loglog
            #ax.set_ylim(1e-11,  1e+3 * 2)  # Test for non-resampling for loglog
        else:
            ax.set_ylim(1e-11,  0.028)
            ax.set_xlim(0.001, 0.6)
    elif("pressure" in param):
        if(plot_loglog):
            ax.set_ylim(1e+1 , 1e+9 * 3) #for loglog
            #ax.set_ylim(1e+1 , 1e+14* 4) # Test for non-ressampling loglog
        else:
            ax.set_ylim(1e+1 , 1e+10 * 0.6) 
            ax.set_xlim(0.001, 0.6)
################




def FFT(data, sampling_rate):
    """
    Perform FFT on velocity data and plot frequency spectrum.
    Parameters:
        velocity_data (array-like): Velocity data for analysis.
        sampling_rate (float): Sampling rate in Hz.
    Returns:
        freqs (numpy.ndarray): Frequencies of the spectrum.
        power (numpy.ndarray): Power of the frequency spectrum.
    """
    #set the sampling_rate to
    sampling_rate = 1.0
    #Detrend the data to remove linear trends
    mean = np.mean(data) 
    data = data - mean
    data = detrend(data)
    #Aplly a window function to reduce spactral leakage
    #window = np.hanning(len(data))
    #data   = data * window
    #Perform FFT
    n    = len(data)
    fft_result = np.fft.fft(data)
    freqs = np.fft.fftfreq(n, d=1/sampling_rate)
    power = np.abs(fft_result)**2

    #Keep only positive frequencies
    positive_freqs = freqs[:n // 2]
    positive_power = power[:n // 2]
    ##sorted find the index of the highest peak
    #p_index = np.argmin(positive_power)
    ##Arrange the data such as that the highest peak corresponds to the maximum energy
    #sorted_index =np.argsort(positive_power)
    #freq_sorted  = positive_freqs[sorted_index]
    #power_sorted  = positive_power[sorted_index]

    #positive_freqs = freq_sorted
    #positive_power = power_sorted
    
    #return value
    return (positive_freqs, positive_power)
    #return (positive_freqs *1000, positive_power)






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
with open("confreq.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
PARAMSDICT     = Fig_params['PARAMSDICT']
################
#Set the SEAFLOOR DEPTH at the location of the ADCP
SEAFLOOR_DEPTH = Fig_params['SEAFLOOR_DEPTH'] # in meter
#Set the ADCP DEPTH
ADCP_DOWN_HEIGTH_FROM_LAKEBED = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']  #in meter
ADCP_UP_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']  #in meter
#Plot Option
Remove_frame   = Fig_params["Remove_frame"]
#Seimic mseed file
mseedFile      =   Fig_params['mseedFile']
mseedFileP      =   Fig_params['mseedFileP']
#Seimic mseed file
Response_File  = Fig_params['Response_File']
#Grab the path of the data
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
#Grab the type of the ADCP data to plot
plot_loglog      = Fig_params['plot_loglog']
#Discharge file name
FILE_DISCH     = Fig_params['FILE_DISCH']
#Grab the Meteo data file
FILE_METEO     = Fig_params['FILE_METEO']
#Grab the temperature from BRB-SOLO-3
FILE_TEMPE     = Fig_params['FILE_TEMPE']
########################################
#Check the user want to make a plot bigger?
plot_loglog    = Fig_params['plot_loglog']
#Grab the starttime
STARTDATE      = Fig_params['STARTDATE']
#Plot Velocity Parameters it True or False
velocity_up    = Fig_params['velocity_up']
velocity_down  = Fig_params['velocity_down']
#Parameter to plot Horizontal velocity Direction (1D) in degrees clockwise from North: degrees_CW_from_N
veldir1D_up    = Fig_params['veldir1D_up']
veldir1D_down  = Fig_params['veldir1D_down']
#Plot the vertical velocity Parameters it True or False
vertical_vel_up = Fig_params['vertical_vel_up']
vertical_vel_down = Fig_params['vertical_vel_down']
#Plot Pressure Parameters it True or False
pressure_up     = Fig_params['pressure_up']
pressure_down   = Fig_params['pressure_down']
###################################
#Plot the depth change of ADCP
depth_adcp_up   = Fig_params['depth_adcp_up']
depth_adcp_down = Fig_params['depth_adcp_down']
###########Surface tracking #########################
surface_tracking_up   = Fig_params['surface_tracking_up']
surface_tracking_down = Fig_params['surface_tracking_down']
##Grab the temperature to visualize
temperature_up   = Fig_params['temperature_up']
temperature_down = Fig_params['temperature_down']
temperature      = Fig_params['temperature']
############################################
##Grab the wind speed paramter to plot
wind_speed       = Fig_params['wind_speed']
##Grab the wind speed direction  paramter to plot
wind_direction   = Fig_params['wind_direction']
#Get the precipation parameter
precipitation   = Fig_params['precipitation']
#Grab the discharge parameter
discharge       = Fig_params['discharge']
#Plot the Average Backscatter?
avg_BS_up       = Fig_params['avg_BS_up']
avg_BS_down     = Fig_params['avg_BS_down']
#Plot the 2D velocity
# Plot the seismic Waveform 
waveform          = Fig_params['waveform']
pressure_waveform = Fig_params['pressure_waveform']
#Grab the bandpass
bandpass          = Fig_params['bandpass']
bandpass_spec     = Fig_params['bandpass_spec']
#get the nbin
n_bin             = Fig_params['n_bin']
#Get the number of range 2D plot 
num_range         = Fig_params['num_range']
#Grab the time list
try:
    date_all_UTC         = pd.date_range(start= STARTDATE, end = STARTDATE, tz = 'UTC')
    #Convert the Switzerland Time into UTC time
except:
     print("Checked the date entring %s"%(STARTDATE))
     exit()
#Create a time list for the started and the
date_list_UTC = [str(it).split()[0] for it in date_all_UTC]
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
        #read ADCP looking downward
        df_down      = dlfn.read(ADCP_FILE_NAME_DOWN) 
        #read ADCP looking upward
        df_up        = dlfn.read(ADCP_FILE_NAME_UP)
        EXT          = os.path.basename(ADCP_FILE_NAME_UP).split("_")[1] 
        #########################################
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
#get the frequency
fs_up        = df_up.fs
fs_down      = df_down.fs
fs_wind      = 1./(10 * 60)
fs_disch     = 1./(10 * 60)
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


#Create the figure
#Grab the figure size
#fig, axs  = plt.subplots(1, nfigs, sharex = False, figsize = fig_size)
if(waveform==True or pressure_waveform==True):
    fig, axs  = plt.subplots(1,nfigs, figsize = fig_size)
else:
    fig, axs  = plt.subplots(nfigs,1, sharex = False, figsize = fig_size)
#
#get the basename
basename_ = os.path.basename(mseedFile) 
basename  = basename_.split(".mseed")[0]
#######################################


####################################################
#################################################################
if(waveform):
    #set the parameter
    param      = "waveform"
    lw         =0.3
    #get the seismogram from mseedFile
    time, data, Title = read_mseed(STARTDATE, mseedFile, Response_File, bandpass,resampling=True)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    fs = 1.0
    #Set the Title
    freqs, power = FFT(data, fs)
    #Plot the spectogram
    Plot_waveforms_freq(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)


if(pressure_waveform):
    #set the parameter
    param      = "pressure_waveform"
    lw         =0.3
    #get the seismogram from mseedFile
    time, data, Title = read_mseed(STARTDATE, mseedFileP, Response_File, bandpass, resampling=True, Pressure = True)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    fs = 1.0
    #Set the Title
    freqs, power = FFT(data, fs)
    #Plot the spectogram
    Plot_waveforms_freq(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)


#Plot the Discharge
if(discharge):
    #set the parameter
    param         = "discharge"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_disch,  date_list_UTC, 'Discharge')
    ####################3
    freqs, power = FFT(data, fs_disch)
    #np.savetxt('Discharge.dat' , data , delimiter="\n", fmt="%.4f")
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

#Doubleched if we need to plot the Velocity
if(velocity_up):
   #set the parameter
   param      = "velocity_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "velocity")
   freqs, power = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, linewidth= 1.2,  **PARAMSDICT)
   #Plot the figure by calling the plotting function
if(velocity_down):
   #set the parameter
   param      = "velocity_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, "velocity")
   freqs, power = FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
   #Plot the figure by calling the plotting function
######################## velocity up ####################
if(vertical_vel_up):
       #set the parameter
       param      = "vertical_vel_up"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg_up , date_list_UTC, param)
       freqs, power = FFT(data, fs_up)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix, _ = PARAMSDICT[param]
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], freqs, data, param, plot_loglog, **PARAMSDICT)

if(vertical_vel_down):
       #set the parameter
       param      = "vertical_vel_down"
       #param      = "vertical_vel"
       #get the vertical velocity ("w") data for the corresponding time
       time, data = Extract_df_list(df_avg_down , date_list_UTC, param)
       freqs, power = FFT(data, fs_down)
       #get the parameter, the index of the corresponding axis, and the color
       _ , color, ix, _ = PARAMSDICT[param]
       data = 1 * data
       #Plot the both velocity on the same figure if there are both set to True by calling the plotting function
       Plot_fig(axs[ix], freqs, power, param, plot_loglog,linewidth=0.8, **PARAMSDICT)
####################################################################
if(veldir1D_up):
   #set the parameter
   param            = "veldir1D_up"
   #get the horizontal velocity data for the corresponding time
   time, data       = Extract_df_list(df_avg_up , date_list_UTC, param)
   freqs, power     = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
   #Plot the figure by calling the plotting function

if(veldir1D_down):
   #set the parameter
   param      = "veldir1D_down"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_down , date_list_UTC, param)
   freqs, power     = FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
   #Plot the figure by calling the plotting function
#set the parameter
if(depth_adcp_up):
   param      = "depth_adcp_up"
   #get the horizontal velocity data for the corresponding time
   time, data = Extract_df_list(df_avg_up , date_list_UTC, "depth")
   freqs, power     = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   data = SEAFLOOR_DEPTH - data
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

if(depth_adcp_down):
   #set the parameter
   param       = "depth_adcp_down"
   #get the horizontal velocity data for the corresponding time
   time, data  = Extract_df_list(df_avg_down , date_list_UTC, "depth")
   freqs, power= FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   data = SEAFLOOR_DEPTH - data
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)


if(surface_tracking_up):
   param        = "surface_tracking_up"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_up , date_list_UTC, "depth")
   freqs, power = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
   #################################################################
if(surface_tracking_down):
   param        = "surface_tracking_down"
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_down , date_list_UTC, "depth")
   freqs, power = FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Get the depth of the ADCP from the surface 
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
   #Plot the figure by calling the plotting function

##Doubleched if we need to plot the Temperature
if(temperature):
   #set the parameter
   param        = "temperature"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_RBR_SOLO_TEMP(data_temp,  date_list_UTC, AVERAGE=True)
   freqs, power = FFT(data, fs_wind)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)


if(temperature_up):
   #set the parameter
   param        = "temperature_up"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_up , date_list_UTC, "temp")
   freqs, power = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

if(temperature_down):
   #set the parameter
   param            = "temperature_down"
   #De-comment this if you want to plot ADCP temparature
   #get the horizontal velocity data for the corresponding time
   time, data   = Extract_df_list(df_avg_down , date_list_UTC, "temp")
   freqs, power = FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
#Plot the Pressure
if(pressure_up):
   #set the pad_echr
   param         = "pressure_up"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg_up , date_list_UTC, "pressure")
   freqs, power = FFT(data, fs_up)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

if(pressure_down):
   #set the pad_echr
   param         = "pressure_down"
   #get the pressure  for the corresponding time
   time, data    = Extract_df_list(df_avg_down , date_list_UTC, "pressure")
   freqs, power  = FFT(data, fs_down)
   #get the parameter, the index of the corresponding axis, and the color
   _ , color, ix, _ = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #convert decibar to bar; 1dbar = 0.1 bar
   data  = 0.1 * data
   Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
##Plot the Turbidity
if(avg_BS_up):
    #set the parameter
    param      = "avg_BS_up"
    #get the Turbidity data for the corresponding time
    time, data = Extract_df_list(df_avg_up,  date_list_UTC, "avg_BS")
    freqs, power  = FFT(data, fs_up)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
##Plot the Turbidity
if(avg_BS_down):
    #set the parameter
    param      = "avg_BS_down"
    #get the Turbidity data for the corresponding time
    time, data = Extract_df_list(df_avg_down,  date_list_UTC, "avg_BS")
    freqs, power  = FFT(data, fs_down)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

if(wind_speed):
    #set the parameter
    param         = "wind_speed"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_UTC, param)
    freqs, power  = FFT(data, fs_wind)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)

if(wind_direction):
    #set the parameter
    param         = "wind_direction"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_UTC, param)
    freqs, power  = FFT(data, fs_wind)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _ = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    Plot_fig(axs[ix], freqs, power, param, plot_loglog, **PARAMSDICT)
#Plot the precipitation
if(precipitation):
    #set the parameter
    param         = "precipitation"
    #get the Dischage data for the corresponding time
    time, data    = Extract_database(d_mteo,  date_list_UTC, param)
    #get the parameter, the index of the corresponding axis, and the color
    _ , color, ix, _= PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    # convert hPa to bar;  1hPa =   0.001 bar
    #Plot_fig(axs[ix], time,  data, param, plot_loglog, **PARAMSDICT)
    Plot_hist(axs[ix], time, data, param, plot_loglog,  **PARAMSDICT)

################################################################################

##Set the ticks on the last x-axis
#axs.set_xlabel('Power (dB/Hz) ', fontsize = 13)
#Set the font size of yticks
#plt.yticks(fontsize=14)
#Set the font size of xticks
#plt.xticks(fontsize=14)
#Aligned all the ylabels of the figure
fig.align_ylabels()
if(pressure_waveform):
    figname   = "Waveform-%s_%s_of_%s_freqs.png"%(EXT, basename, STARTDATE)
    fig.subplots_adjust(wspace=0.3)
else:
    figname   = "ADCP-%s_%s_of_%s_freqs.png"%(EXT, basename, STARTDATE)
    axs[nfigs -1].set_xlabel('Frequency (Hz) ', fontsize = 13)
##Space between the subplots
if(nfigs <= 5):
    plt.subplots_adjust(hspace = fig_space)
    fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
elif(nfigs > 5):
    fig_space = fig_space + (1.0/100) * nfigs
    plt.subplots_adjust(hspace = fig_space)
    #Aligned all the ylabels of the figure
    fig.align_ylabels()
    #Save the figure
    fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
