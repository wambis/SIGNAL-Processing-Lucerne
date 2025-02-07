import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates
import re, pandas as pd 
import dolfyn as dlfn
import datetime
import json, yaml
from yaml.loader import SafeLoader
from obspy import read, read_inventory
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
#################
from scipy.fft import fft, fftfreq
#from scipy.signal import spectrogram
import scipy
from scipy import interpolate
#################################################
from matplotlib.ticker import FormatStrFormatter
#
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

def inter1D(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

def inter1DD(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)

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

def check_if_string(var):
    return isinstance(var, str)

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





#def Plot3DBackscatter(ax, df,TARGETDATE, height_adcp, depth_idx, ADCP_DOWN=False, **PARAMSDICT):
def Plot3DBackscatter(ax, df,TARGETDATE, height_adcp, depth_idx, ADCP_DOWN=False, VERTICAL_SLIDES =None):
    #get the backscatter for the target date
    backscat     = df.amp.sel(time = TARGETDATE)
    #get the paramter
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = backscat.range.data
    time         = backscat.time.data
    beam_indx    = backscat.beam.data
    #Check the if you're plotting the upward, downlard
    if(ADCP_DOWN):
        depths   = abs( (height_adcp + blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
    else:
        depths   = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
    #Dimensions: time x depth x beam
    n_time  = len(time)
    n_depth = len(depths)
    n_beam  = len(beam_indx)
    ################################# 
    beams   = beam_indx  # Beam index
    #get the time in the matplotlib type 
    times = mdates.date2num(time)
    # Highlighted slice (e.g., at depth index 10)
    highlight_depth = depths[depth_idx]
    #Plot the verticals slice
    if(VERTICAL_SLIDES):
        T, B = np.meshgrid(beams, times, indexing='ij')
        D    = np.full_like(T, highlight_depth)
        #loop over the depths
        for d in range(n_depth):
            slice_data  = backscat[:, d, :]
            depth_layer = np.full_like(T, depths[d])
            #ax.plot_surface(T, B, depth_layer, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
            ax.plot_surface(B,T, depth_layer, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
                            shade=False, alpha=0.7, edgecolor='none')
        
        #Highlighted slice
        highlight_data = backscat[:, depth_idx, :]
        #ax.plot_surface(T, B, D, rstride=10, cstride=1, facecolors=cm.PuBu(highlight_data / highlight_data.max()),
        ax.plot_surface(B, T, D, rstride=10, cstride=1, facecolors=cm.PuBu(highlight_data / highlight_data.max()),
                        shade=False, alpha=0.9, edgecolor='k')
        
    #Plot the horizontal slice
    else:
            T, B = np.meshgrid(beams, times, indexing='ij')
            D    = np.full_like(T, highlight_depth)
            #Plot 3D backscatter slices at each depth
            for d in range(n_depth):
                slice_data = backscat[:, d, :]
                depth_layer = np.full_like(T, depths[d])
                ax.plot_surface(T, depth_layer, B, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
                                shade=False, alpha=0.7, edgecolor='none')
            #Highlighted slice
            highlight_data = backscat[:, depth_idx, :]
            ax.plot_surface(T, D, B, rstride=10, cstride=1, facecolors=cm.plasma(highlight_data / highlight_data.max()),
                            shade=False, alpha=0.9, edgecolor='k')




#Open the configuration file
with open("confbck.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)

#Grab the path of the data
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
#Discharge file name
#Grab the Meteo data file
#############################################
#Check the user want to make a plot bigger?
plot_bigger    = Fig_params['plot_bigger']
#Grab the starttime
TARGETDATE     = Fig_params['TARGETDATE']
#Grab the dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#get the type of slice
VERTICAL_SLIDES  = Fig_params['VERTICAL_SLIDES']
##########################################
ADCP_DOWN      = Fig_params['ADCP_DOWN']
######## Seimic Data file########################
mseedFile      = Fig_params['mseedFile']
xmlFile        = Fig_params['xmlFile']
bandpass       = Fig_params['bandpass']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED     = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']

#get the depth index the user want to plot
depth_idx    = Fig_params['depth_idx'] 


if(ADCP_DOWN):
    df           = dlfn.read(ADCP_FILE_NAME_DOWN) 
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    depths       = abs( (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
elif(ADCP_UP):
    df           = dlfn.read(ADCP_FILE_NAME_UP) 
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    depths       = (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
else:
    print("print provide the ADCP file %s "%(ADCP_FILE_NAME_DOWN))

    
#get the backscatter
param = "Backscatter1D"
tad  , d_amp = Extract_ADCPData(df,  TARGETDATE, param)
# Convert pandas datetime to matplotlib date format
tadd    = mdates.date2num(tad)

#get the seismogram from mseedFile
time, tr, Title   = read_mseed(TARGETDATE, mseedFile, xmlFile, bandpass)
##
ts, Power         = seismic_FFT(tr, TARGETDATE)
#Make the interpolation of the seismic energy to match the backscatter
tsnew, Powernew  =  inter1DD(ts,  Power, tadd)



#Create 3D plot
#fig = plt.figure(figsize=(12, 14))
fig = plt.figure(figsize=(8, 10))
# Set height ratios to 1:1 for equal sizes
#gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])  
#gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1],  width_ratios=[1], wspace=0.05, hspace=0.05)  
gs_main = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 0, 0])  # Third and fourth row empty
#gs = gridspec.GridSpec(2, 1, width_ratios=(2, 1), height_ratios=(1, 2), wspace=0.05, hspace=0.05)  
#gs = gridspec.GridSpec(2, 1, width_ratios=[1], wspace=0.05, hspace=0.05)  
#Add for  3D plot
#ax  = fig.add_subplot(211, projection='3d')  # Top
ax  = fig.add_subplot(gs_main[0], projection='3d')  # Top
#Create a nested GridSpec to reduce 1D plot width
gs_sub = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = gs_main[1], width_ratios=[1, 4, 1])
#gs_sub = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = gs_main[1], width_ratios=[1., 4, 0.5])
#gs_sub = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = gs_main[1], width_ratios=[0.5, 4, 1])
#Add axis for 1D plot
#ax2 = fig.add_subplot(212)  # Bottom
ax2 = fig.add_subplot(gs_sub[1], sharex=ax)  # Bottom
#Set aspect ratio to be the same
#ax2.set_aspect(aspect='auto')
#ax.set_position([1, 1, 1, 1])  # [left, bottom, width, height]
#ax2.set_position([0.5, 0.5, 0.1, 0.1])  # [left, bottom, width, height]

# Set the aspect ratio of the 3D plot to match the 1D plot
#ax.set_box_aspect([1, 1, 1])
#ax2.set_box_aspect([1, 1, 0.7])
#Labels and view adjustment
fsize= 12
lpad = 12
lsize = 11
pad = 10
#set this for time limit
tstart = datetime.datetime(2023,10,23,0,0)
tend   = datetime.datetime(2023,10,24,0,0)
#Check if the user want to plot the vertical slice
if(VERTICAL_SLIDES):
    #figname = "Backscatter_3D_Vertical_Slices.png"
    figname = "Backscatter_3D_Vertical_Slices.pdf"
    if(ADCP_DOWN):
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_DOWN_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    else:
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_UP_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_ylabel("Beams", fontsize=fsize, rotation=100,labelpad =lpad)
    ax.set_ylabel("Beams", fontsize=fsize)
    ax.set_zlabel("Height above lakebed (m)", fontsize=fsize, labelpad =lpad)
    ax.set_zlim(top=max(depths))
    #ax.set_zlim(min(depths), max(depths))
    #number of vertical lines for grid
    locations    = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #control the visibility of the grid
    ax.set_box_aspect([1,1,1])
    ax.grid(False)
    ax.set_facecolor((1,1,1,0))
    ax.xaxis.pane.fill = False  #Remove pane background
    ax.yaxis.pane.fill = False  #Remove pane background
    ax.zaxis.pane.fill = False  #Remove pane background
    ########################################3
    #Set ylim of x-axis
    ax.set_xlim([tstart, tend])
    #ax.set_ylabel("Time (hour:minute)",fontsize=fsize, labelpad =lpad)
#    ax.set_xlabel("Time (hour:minute)",fontsize=10, labelpad =12)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ##################################################
    #ax.tick_params(axis ='y', labelsize= lsize, pad=1.5)
    ax.tick_params(axis ='y', labelsize= lsize, pad=1.5)
    ax.tick_params(axis ='x', labelsize= 9, labelrotation=310,zorder=12, pad=0.5)
    ax.tick_params(axis ='z', labelsize= lsize,zorder=12)
    #ax.set_title("3D ADCP Backscatter with Highlighted Depth Slice")
    #ax.view_init(elev=20, azim=-30)  # Adjust view angle
    ax.view_init(elev=9.5, azim=-82)  # Adjust view angle
else:
    figname = "Backscatter_3D_Horizontal_Slices.png"
    pad_ticks_h  = 1.2
    if(ADCP_DOWN):
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_DOWN_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    else:
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_UP_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_xlabel("Beam index", fontsize=fsize, labelpad =lpad)
    ax.set_xlabel("Beams", fontsize=fsize, labelpad =lpad)
    ax.set_ylabel("Height above lakebed (m)", fontsize=fsize, labelpad =lpad)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.zaxis.set_major_locator(locations)
    #########################################3
    #Set zlim of x-axis
    ax.set_zlim([tstart, tend])
    ax.set_zlabel("Time (hour:minute)",fontsize=fsize, labelpad =lpad)
    ax.zaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ##################################################
    ax.tick_params(axis ='x', labelsize= lsize, pad= pad_ticks_h)
    ax.tick_params(axis ='y', labelsize= lsize, pad= pad_ticks_h)
    ax.tick_params(axis ='z', labelsize= lsize)
    #ax.set_title("3D ADCP Backscatter with Highlighted Depth Slice", loc="center", pad = -15.0)
    #ax.view_init(elev=30, azim=-45)  # Adjust view angle
    ####################
    #ax.view_init(elev=20, azim=-30)  # Adjust view angle
    #ax.view_init(elev=20, azim=10)  # Adjust view angle
    ax.view_init(elev=-100, azim=270)  # Adjust view angle

# Create a second y-axis sharing the same x-axis
ax3 = ax2.twinx()
#plot the 1D plot here on the ax2 axis
fig.canvas.draw()  # Ensure the figure is fully rendered before getting positions
# Adjust position of ax2 to shift it left
pos3D = ax.get_position()  # Get bounding box of 3D plot
pos1D = ax2.get_position()  # Get bounding box of 1D plot

#new_pos = [pos3D.x0 + 0.06, pos1D.y0 + 0.05, pos1D.width, pos1D.height]  # Align left with 3D plot
new_pos = [pos3D.x0,  pos1D.y0 , pos1D.width, pos1D.height]  # Align left with 3D plot
#new_pos = [pos3D.x0 - pos3D.x0 *0.9 ,  pos1D.y0 , pos1D.width , pos1D.height]  # Align left with 3D plot
new_3D = [pos3D.x0 +2,  pos1D.y0 , pos1D.width , pos1D.height]  # Align left with 3D plot
ax.set_position(new_3D)  # Apply new position
ax2.set_position(new_pos)  # Apply new position
ax3.set_position(new_pos)  # Apply new position



#plot the Seismic energy
ax2.plot(tadd, Powernew, color='k', lw=2, label='Seismic Energy')
ax2.plot(tadd, Powernew, color='k', lw= 14, alpha =0.3 )

#Plot the back scatter
bs_color = '#0a481e'
#ax3.plot(tadd, d_amp, color='r', label='Backscatter')
ax3.plot(tadd, d_amp, color=bs_color, lw =2., label='Echo')
ax3.plot(tadd, d_amp, color=bs_color, lw=14, alpha = 0.3)

########################################
locations    = verical_grid()
ax2.xaxis.set_major_locator(locations)
ax2.tick_params(axis ='x', labelsize= 9, labelrotation=310,zorder=12, pad=0.5)
#set the labels
ax2.set_ylabel("Seismic energy (dB)", color ='k', labelpad = pad, fontsize=fsize)
ax3.set_ylabel('Echo (counts)', color=bs_color, labelpad = pad, fontsize=fsize)

ax2.set_xlim([tstart, tend])
ax2.set_xlabel("Time (hour:minute) on %s"%(TARGETDATE),fontsize=fsize, labelpad =pad)
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))



#Add the colorbar and save the figure
plt.colorbar(cm.ScalarMappable(cmap='Spectral'), ax=ax, label='Echo Intensity', pad =0.04, extend ='both', shrink=0.5)
#ax.colorbar(cm.ScalarMappable(cmap='Spectral'), ax=ax, label='Backscatter Intensity', pad =0.04, extend ='both', shrink=0.5)
# Adjust layout to prevent overlap
#plt.tight_layout()  
# Improve layout to ensure same size
#fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.2)
plt.subplots_adjust(hspace=0.1)  # Reduce vertical spacing globally
plt.savefig(figname, bbox_inches = 'tight')
