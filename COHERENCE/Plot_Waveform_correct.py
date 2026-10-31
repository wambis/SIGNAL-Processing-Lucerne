##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
from obspy import read, read_inventory, Stream
import numpy as np
from scipy.signal import welch
from obspy import read
from scipy import signal
import os, re
import pandas as pd
import matplotlib.dates as mdates
import datetime

#Set parameters
#xy_coords   = (0.85, 0.95)
xy_coords   = (0.85, 0.97)
#xy_coords_D = (0.55, 0.85)
#set the fontsize
Figsize= 14
pad    = 10
#define text size
f_size = 14
####for ticks #####
#fticks = 14
fticks = 18


##################################################################
#import numpy as np
#from scipy import signal


#Read the mseed file
def read_mseed(FILE, xmlFile, bandpass, Pressure = None):
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

    #####################################
    st.remove_sensitivity(inventory=inv)
    st.filter("bandpass", freqmin=f2, freqmax=f3)

    #if(Pressure):
    #    st.remove_response(inventory=inv, pre_filt=pre_filt, output="DEF",water_level=60)
    #else:
    #    #Remove instrument response
    #    st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=30)
    #    #st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
    #    #st.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",water_level=60)

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
        if(component=="1"):
            component ="Z"
        elif(component=="2"):
            component ="HH2"
        elif(component=="3"):
            component ="HH2"
        #Title         = "%s %s %s %s  %s  %s"%(network, station, channel, dT, stime.split(".")[0], endtime.split(".")[0])
        Title         = "%s.%s.%s"%(network, station, component)
        ##Resampling the data
        tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
        #################################
        t0 = tr.stats.starttime.datetime        # Python datetime
        #sampling_rate = tr.stats.sampling_rate
        data          = tr.data
        #time          = tr.times()
        Date_of_Day   = str(tr.stats.starttime.date)
        #The sampling rate is 1 second
        date_rng      = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
        #change the date_rng to numpy
        #change the date_rng to numpy
        time         = date_rng.to_numpy()
    return(t0, time, data, Title)




def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations


def remove_pressure_coherence(Z, P, fs, nperseg=4096):
    """
    Remove the part of the vertical seismic signal Z that is coherent with pressure P.
    Computes coherence using scipy.signal.coherence inside the function.

    Parameters
    ----------
    Z : array
        Vertical seismic component
    P : array
        Pressure signal
    fs : float
        Sampling frequency (Hz)
    nperseg : int
        Segment length for spectral estimation

    Returns
    -------
    Z_clean : array
        Vertical component with pressure-coherent energy removed
    Z_pred : array
        Pressure-induced prediction that was subtracted
    H_interp : array
        Transfer function used for subtraction
    freqs : array
        Frequency vector for FFT
    f_coh : array
        Frequencies for coherence
    Cxy : array
        Coherence between Z and P
    """
    #--- 1. Compute coherence (your requirement) ---
    f_coh, Cxy = signal.coherence(Z, P, fs=fs, nperseg=nperseg)

    #--- 2. Compute cross-spectral and auto-spectral densities ---
    f, P_ZP = signal.csd(Z, P, fs=fs, nperseg=nperseg)
    f, P_PP = signal.welch(P, fs=fs, nperseg=nperseg)

    #--- 3. Wiener filter transfer function ---
    H       = P_ZP / P_PP

    #--- 4. FFT of pressure signal ---
    N       = len(Z)
    freqs   = np.fft.rfftfreq(N, 1/fs)
    P_fft   = np.fft.rfft(P, n=N)

    #--- 5. Interpolate transfer function to FFT frequencies ---
    H_interp  = np.interp(freqs, f, H.real) + 1j * np.interp(freqs, f, H.imag)

    #--- 6. Predict pressure-induced vertical motion ---
    Z_pred_fft= H_interp * P_fft
    Z_pred    = np.fft.irfft(Z_pred_fft, n=N)

    #--- 7. Subtract to obtain cleaned vertical component ---
    Z_clean   = Z - Z_pred

    return Z_clean, Z_pred, H_interp, freqs, f_coh, Cxy

#def start_endtime(time_serie, ListDates):
#    #print(time_serie)
#    #exit()
#    #minimum time
#    tstart     = time_serie.min()
#    #tstart     = time_serie[0]
#    #print(tstart)
#    #exit()
#    #Add one day to the minimum time
#    #tend      = tstart + pd.to_timedelta(len(ListDates), unit= 'D')
#    tstart = pd.Timestamp(tstart)
#    tend   = tstart + pd.to_timedelta(len(ListDates), unit='D')
#
#    #tend       = time_serie[-1]
#    #Convert to minimum
#    tend       = tend.to_numpy()
#    return (tstart, tend)

def start_endtime(time_serie, ListDates):
    #Convert minimum time to pandas Timestamp
    tstart = pd.to_datetime(time_serie.min())

    # Add number of days in ListDates
    tend = tstart + pd.to_timedelta(len(ListDates), unit='D')

    # Convert to numpy datetime64 if needed
    tend = tend.to_numpy()

    return (tstart, tend)



################################################
def plot_raw_vs_corrected(Z, Z_clean, t0, fs, Title_z):
    """
    Plot the raw and corrected vertical seismic signals on the same graph.

    Parameters
    ----------
    Z : array
        Raw vertical seismic signal
    Z_clean : array
        Corrected vertical seismic signal
    fs : float
        Sampling frequency (Hz)
    """

    N = len(Z)
    t = np.arange(N) / fs
    # Convert each time value (seconds) into a datetime
    t_datetime = [t0 + datetime.timedelta(seconds=float(tt)) for tt in t]
    # convert to matplotlib time
    #t = mdates.date2num(t0 + t)                
    # Convert to matplotlib time numbers
    t = mdates.date2num(t_datetime)
    #############################
    #get the start and the endtime
    ListDates = [str(t0.date())]
    tstart, tend =  start_endtime(t, ListDates)
    ###########################################
    ylabel = 'ACC (µm/s²)'
    nfigs  = 2
    lw     = 0.25
    c      = 'k'
    ######################
    fig_size= (12, 8)
    ##################Create a figure#######
    fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)

    #plt.figure(figsize=(12, 5))
    axs[0].plot(t, Z * 1e+6,       lw=lw, c=c,  label="Raw signal", alpha=0.7)
    axs[1].plot(t, Z_clean * 1e+6, lw=lw, c=c, label="Corrected signal", alpha=0.7)
    ##########################

    #Set the tile
    axs[0].set_title(
    Title_z     ,
    fontsize=fticks,
    color="black",
    loc="center"   # 'left', 'center', 'right'
    )

    #Add the legend
    locations = verical_grid()
    for i in range(len(axs)):
        axs[i].xaxis.set_major_locator(locations)
        axs[i].grid(visible = True, axis = "both", alpha = 0.7)
        axs[i].legend(loc="upper left")
        #Customize y-axis tick label font size
        #axs[i].tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
        axs[i].tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
        #ax.legend(loc="upper left", fontsize='large',  bbox_to_anchor=(0, 1), frameon=True, shadow=False)
        #set limit
#        axs[i].set_ylim(-0.3, 0.3) #Good for Discharge
        #costumize the y-axis
        axs[i].set_ylabel(ylabel, labelpad = pad, fontsize = fticks)
        #axs[i].set_xlim(tstart, tend)

    
    #Set ylim of x-axis
    #disable axis
    axs[0].set_xticklabels([])

    #axs[0].set_ylabel(ylabel, labelpad = pad, fontsize = fticks)
    #axs[1].set_ylabel(ylabel, labelpad = pad, fontsize = fticks)
    ##############################################
    #axs[1].set_xlabel('Time (s)', labelpad = pad, fontsize = fticks)
    axs[-1].set_xlabel('Time (hour:minute) on %s'%(ListDates[0]), labelpad =12, fontsize = fticks)
    #save the figure
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    figname = "waveforms-correct-%s.pdf"%(ListDates[0])

    fig.subplots_adjust(hspace=0.07)
    #Save the figure
    plt.savefig(figname, bbox_inches = 'tight', dpi = 300)




#
#---------------------------------------------------------


#Example: Reading MSEED files
#north    = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed'
#east     = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed'
#Wind event
FileH1 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH1-2024-03-29.mseed'
FileH2 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH2-2024-03-29.mseed'
FileH3 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH3-2024-03-29.mseed'
FileP  = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HDH-2024-03-29.mseed'


#FileH1 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH1-2024-01-13.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH2-2024-01-13.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH3-2024-01-13.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HDH-2024-01-13.mseed'

###Discharge event
#FileH1 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH1-2023-09-22.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH2-2023-09-22.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH3-2023-09-22.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HDH-2023-09-22.mseed'

#get the response file
Response_File= 'RESPOSE-FILE.XML/xj-fdsn.xml'

#Just removed instrument response
#bandpass= '8-10-80-90' #Good for correction of Discharge day 2023-09-22 for JGR
bandpass= '5-10-90-100' #Good for correction of wind day  2024-03-29 for JGR
# ---------------------------------------------------------
# 1. Load your three-component seismic data
# ---------------------------------------------------------

# Replace with your own file or stream

#Sampling frequency
#fs = trZ.stats.sampling_rate
fs=1.0

t0, time, Z, Title_z =read_mseed(FileH1, Response_File, bandpass)
t0, time, P, Title_p =read_mseed(FileP,  Response_File, bandpass, Pressure=True)

#Z= trZ.data
#P= trP.data
#Compute the corrected signal and the coherence 
Z_clean, Z_pred, H, freqs, f_coh, Cxy = remove_pressure_coherence(Z, P, fs)
###################################
plot_raw_vs_corrected(Z, Z_clean,t0, fs, Title_z)
#####################################################
