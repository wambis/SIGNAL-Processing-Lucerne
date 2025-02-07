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
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
from yaml.loader import SafeLoader
from obspy import read, read_inventory
import matplotlib.ticker as mticker
#################
from scipy.fft import fft, fftfreq
#from scipy.signal import spectrogram
import scipy
from scipy import interpolate


######Oceanic colors
import cmocean

#paramerters
AphaDict = {0:'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            5: 'f', 6: 'f', 7:'h', 8:'g', 9: 'j',10:'k',
            11:'l', 12: 'm',13: 'n',14:'o', 15:'p', 16: 'q'}

##############################################


Figsize= 14
pad    =10

#Define the grid
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

def inter1DD(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)


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
    elif(param=="Perturbation"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
        Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist])
        T           = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
        DList = np.copy(Matrix2D)
        #for i in range(Matrix2D.shape[0]):
        for i in range(len(DList)):
            #print(Matrix3D[beam_indx][i])
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


###############################################
def Extract_BS2D(df, Tlist,  param, nrange=None):
        if(nrange== None):
            r        = df.range.data
        else:
            r        = df.range.data[:nrange]
        if(param=="backscatter"):
                #make the average on all the 4 beams, NB the key word here is amp
                #Get the range, altitude
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(df.amp, axis = 0)
                #Loop over the list of the time
                #mean_val = np.mean(np.mean(df_beam_avg, axis = 1), axis = 0)
                #Loop over the list of the time
                #Matrix3D    = np.concatenate([df.amp.sel(time = ti, range = r) for ti in Tlist], axis =1 )
                Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist])
                DList = np.copy(Matrix2D)
                #for i in range(Matrix2D.shape[0]):
                for i in range(len(DList)):
                        #print(Matrix3D[beam_indx][i])
                        mean_val = np.mean(DList[i])
                        DList[i] = [(item -mean_val)/mean_val for item in DList[i]]
                        #List_l = [(item -mean_val) for item in Matrix3D[beam_indx][i]]
                        #DList.append(List_l)
                #transform DList into an array 
                #DList       = np.asarray(DList)[::-1]
                DList       = np.asarray(DList)
                #print(DList)
                #print(DList.shape)
                #exit()
                time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
                #return (time, r, Matrix3D[beam_indx])
                return (time, r, DList)


def Extract_RBR_SOLO_TEMPERATURE(DataDict,  Tlist, AVERAGE=False):
    #The data should be in the dictionary, where the keys are the depths
    DATA2D     = []
    DEPTHS     = []
    #get the time you need to plot
    TIMES      = set([time for depth in DataDict.keys() for time in DataDict[depth] if(time.split()[0] in Tlist)])
    #####################################
    TIMES      = sorted(np.asarray(list(TIMES)))
    #loop over the depths and the time to form a matrix
    #depthss = [d for d in DataDict.keys()]
    #depths_up   = [70, 60, 50, 37.5]
    #depths_down = [20, 17.5, 15, 5 ]

    depths_up   = [40,  55]
    depths_down = [5, 12.5]
    ##################################
#    depths_up   = [70]
#    depths_down = [20]
    #['70', '67.5', '65', '62.5', '60', '57.5', '55', '52.5', '50', '47.5', '45', '42.5', '40', '37.5', '35', '32.5', '30', '27.5', '25', '22.8', '20', '17.5',  
    #'15', '12.5', '10', '7.5', '5']
    #print(depthss)
    DICT_UP   = dict()
    DICT_DOWN = dict()
    #exit()
    #for  depth in DataDict.keys():
    for  depth in depths_up:
        #get the temperature at a specific depth
        temp   = [DataDict[str(depth)][time] for time in TIMES]
        #####################################
        DICT_UP[depth] = (pd.to_datetime(TIMES), temp)
    ########### loop over the depth down
    for  depth in depths_down:
        #get the temperature at a specific depth
        temp   = [DataDict[str(depth)][time] for time in TIMES]
        #####################################
        DICT_DOWN[depth] =(pd.to_datetime(TIMES), temp)
    #Trans formed the List into# 2D array
    return(DICT_UP, DICT_DOWN)




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










def Plot_fig2D(ax, ix, fig_object, df, time, data2D, height_adcp,  ADCP_DOWN=False, rframe= None):
    #speed up figure plot
    plt.style.use('fast')
    data2D = data2D * 100
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    ylabel = "dln(Echo)";  
    #color = "jet" 
    #color = cmocean.cm.thermal
    #color = cmocean.cm.haline
    #color = cmocean.cm.deep_r
    #color = cmocean.cm.speed_r
    #color = cmocean.cm.rain_r
    #color = cmocean.cm.curl
    color = "RdYlBu_r" 
    ############################
    # Create custom colormap
    #color = LinearSegmentedColormap.from_list("custom_blue_yellow", colors, N=256)
    #ylabel = 'd'+r'$\ln(Echo)$';  
    ylabel = 'd'+r'$\ln$'+'(Echo)'  
    #color = "jet" 
    #ylabel = "Backscatter";  color = "jet" 
    vmind         =min(np.min(data2D, axis = 1))
    vmaxd         = max(np.max(data2D, axis = 1))
    ########################
    #vmind = -50; vmaxd =+50
    vmind = -40; vmaxd =+50
    #vmind = -40; vmaxd =+40
#    vmind = -50; vmaxd =+70
    #vmind = -70; vmaxd =+60
    #vmind = -0.5; vmaxd =+0.5
    #vmind = -0.8; vmaxd =+0.8
    if("velocity" in ylabel):
        vmind = 0.0; vmaxd = 0.4
    #####Test sur L'ADCP orientation
    if(ADCP_DOWN):
        #invert the y-axis
        #ax.invert_yaxis()
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
    tstart, tend  =  start_endtime(time)
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
    if(ADCP_DOWN):
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
    ax.set_xlim(tstart, tend)
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])




#Open the configuration file
with open("conf2Dbck.yaml") as Fym:
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
#get the date
STARTDATE    = Fig_params['STARTDATE']
ENDDATE      = Fig_params['ENDDATE']
#Grab the dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#get the type of slice
VERTICAL_SLIDES  = Fig_params['VERTICAL_SLIDES']
##########################################
ADCP_DOWN      = Fig_params['ADCP_DOWN']
ADCP_UP        = Fig_params['ADCP_UP']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED     = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']



######## Seimic Data file########################
mseedFile      = Fig_params['mseedFile']
xmlFile        = Fig_params['xmlFile']
bandpass       = Fig_params['bandpass']
###############################################


try:
    df_down = dlfn.read(ADCP_FILE_NAME_DOWN)
    df_up   = dlfn.read(ADCP_FILE_NAME_UP)
except:
    print("print provide the ADCP files:  %s and %s  "%(ADCP_FILE_NAME_UP, ADCP_FILE_NAME_DOWN))


try:
    date_all = pd.date_range(start= STARTDATE, end = ENDDATE)
except:
     print("Checked the date entring %s  %s"%(STARTDATE, ENDDATE))
     exit()
#Create a time list for the started and the
date_list = [str(it).split()[0] for it in date_all]



#get the backscatter
#param = "Backscatter1D"
param = "Perturbation"
tad  , d_amp = Extract_ADCPData(df_down,  STARTDATE, param)
#tad  , d_amp = Extract_ADCPData(df_up,  STARTDATE, param)
d_amp  = d_amp * 100
# Convert pandas datetime to matplotlib date format
tadd    = mdates.date2num(tad)

#get the seismogram from mseedFile
time, tr, Title   = read_mseed(STARTDATE, mseedFile, xmlFile, bandpass)
##
ts, Power         = seismic_FFT(tr, STARTDATE)
#Make the interpolation of the seismic energy to match the backscatter
tsnew, Powernew  =  inter1DD(ts,  Power, tadd)

#tstart = datetime.datetime(2023,10,23,0,0)
#tend   = datetime.datetime(2023,10,24,0,0)
tstart, tend = start_endtime(time)




#get the backscatter
if(ADCP_UP==True and ADCP_DOWN==True):
    nfigs = 3
elif(ADCP_UP==True and ADCP_DOWN==False or ADCP_UP==False and ADCP_DOWN==True):
    nfigs = 1

#grab the temperature file
FILE_TEMPE     = Fig_params['FILE_TEMPE']

#Grab the space between the figure
fig_space = Fig_params['fig_space']

#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

#Create the figure
fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)


#Grab 2D Temperature
#get the temperature data from RBR-SOLO-3
# Opening JSON file
ftmp   = open(FILE_TEMPE)
# returns JSON object as
# a dictionary
data_temp = json.load(ftmp)
#get the horizontal velocity data for the corresponding time
fsize= 12
lpad = 12
lsize = 11
pad = 10
#set this for time limit
#######################################
if(ADCP_UP==True and ADCP_DOWN==True):
        param = 'backscatter'
        #for bindx in range(4):
        time, rup, Matrix2D_up     = Extract_BS2D(df_up,  date_list,  param, nrange=None)
        time, rdown, Matrix2D_down = Extract_BS2D(df_down,  date_list,  param, nrange=None)
        #plot the 2D Backscatter
        #plot the upward looking ADCP
        Plot_fig2D(axs[0], 0, fig, df_up, time, Matrix2D_up, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, temps_up, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, DICT_UP, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #plot the downward looking ADCP
        Plot_fig2D(axs[1], 1, fig, df_down, time, Matrix2D_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        #Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, temps_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        #Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, DICT_DOWN, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        # Optional: Set the same limits or aspect ratio
        #plot the Seismic energy
        #Create a second y-axis sharing the same x-axis
        ax2     = axs[2]
        ax3     = ax2.twinx()
        ##############
        ax2.plot(tadd, Powernew, color='k', lw=2, label='Seismic Energy')
        ax2.plot(tadd, Powernew, color='k', lw= 14, alpha =0.3 )
        
        #Plot the back scatter
        #bs_color = '#0a481e'
        bs_color = 'r'
        #bs_color = '#0B2E5B'
        ylabel_echo = 'd'+r'$\ln$'+'(Echo)'
        #ylabel_echo = 'Echo (counts)'
        #ax3.plot(tadd, d_amp, color='r', label='Backscatter')
        ax3.plot(tadd, d_amp, color=bs_color, lw =3., label='Echo')
        ax3.plot(tadd, d_amp, color=bs_color, lw=14, alpha = 0.3)
        
        ########################################
        locations    = verical_grid()
        ax2.xaxis.set_major_locator(locations)
        #ANotate
        ax2.annotate(
            AphaDict[2],     # Text for the annotation
            xy=(0.0089, 0.93),                   # Coordinates (relative to axes) for the annotation
            #xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
            xycoords='axes fraction',    # Use axes fraction for positioning
            fontsize=16,                 # Font size for the annotation
            #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
            bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
        )
        #Make the grid of the axis
        ax2.grid(visible = True, axis = "both", alpha = 0.4)
        #ax2.tick_params(axis ='x', labelsize= 9, labelrotation=310,zorder=12, pad=0.5)
        ax2.tick_params(axis='both', which='major', labelsize=Figsize)
        ax3.tick_params(axis='y', which='major', labelsize=Figsize)
        #set the labels
        ax2.set_ylabel("Seismic energy (dB)", color ='k', labelpad = pad, fontsize=fsize)
        #ax3.set_ylabel('Echo (counts)', color=bs_color, labelpad = pad, fontsize=fsize)
        ax3.set_ylabel(ylabel_echo, color=bs_color, labelpad = pad, fontsize=fsize)
        ######background color #######
        #ax3.set_facecolor(bs_color)
        #remove label
        #ax3.set_yticks([])
        ###set limit
        ax2.set_xlim([tstart, tend])

elif(ADCP_UP==True and ADCP_DOWN==False):
        param = 'backscatter'
        time, r, Matrix2D  = Extract_BS2D(df_up,  date_list,  param)
        Plot_fig2D(axs, fig, df_up, time, Matrix2D, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)

elif(ADCP_UP==False and ADCP_DOWN==True):
        param   = 'backscatter'
        time, r, Matrix2D   = Extract_BS2D(df_down, date_list,  param)
        #plot the downward looking ADCP
        Plot_fig2D(axs, fig, df_down, time, Matrix2D, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        print(Matrix2D.shape)
#################### start ###########################
plt.subplots_adjust(hspace = fig_space)


#######
if(nfigs>1):
    #axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #Write on the ax2-axis 
    axs[2].set_xlabel("Time (hour:minute) on %s"%(STARTDATE),fontsize=fsize, labelpad =pad)
else:
    axs.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

#figname = "Figure_Backscatter2D_Seismic_Energy_%s.png"%(STARTDATE)
figname = "Figure_Backscatter2D_Seimic_Energy_%s.pdf"%(STARTDATE)
#Save the figure
plt.savefig(figname, bbox_inches = 'tight')
