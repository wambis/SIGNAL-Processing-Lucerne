##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
from scipy.signal import welch
from obspy import read, read_inventory
from obspy.signal import PPSD
from scipy import signal
from pathlib import Path as Pathdir
import json, yaml
from yaml.loader import SafeLoader
import os,re, sys, math
from obspy import Trace, Stream
##############################
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
########################
from scipy.signal import savgol_filter


########################################################
AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}
#size for anotation
#xy_coords = (0.71, 0.95)
xy_coords   = (0.6, 0.95)
xy_coords_D = (0.6, 0.86)
#set the fontsize
Figsize= 14
pad    = 20
#define text size
f_size = 14
####for ticks #####
fticks = 14





#Example: Reading MSEED files
#north    = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed'
#east     = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed'

#FileH1 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH1-2024-03-29.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH2-2024-03-29.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH3-2024-03-29.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HDH-2024-03-29.mseed'




#def read_and_merge_trace(file_list):
#    #set the chunk size
#    chunk_size = 20
#    merged     = []
#    for i in range(0, len(file_list), chunk_size):
#        chunk  = file_list[i:i+chunk_size]
#        ####################################### 
#        with ProcessPoolExecutor() as exe:
#            arrays = list(exe.map(load_data, chunk))
#        #add into the list 
#        merged.append(np.concatenate(arrays))
#    #return the merged traces 
#    merged  = np.concatenate(merged)
#    return merged
#

def load_data(file):
    return read(file)[0]


def find_resonance_frequency(f, hvsr):

    idx = np.argmax(hvsr)
    f0 = f[idx]
    amp = hvsr[idx]

    return f0, amp

#def calculate_hvr_psd(file_Z, file_1, file_2, xmlfile, periode=None):
#    #read the mseedfile
#    stZ= read(file_Z)
#    st1= read(file_1)[0]
#    st2= read(file_2)[0]
#    data1  = np.nan_to_num(st1.data.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
#    data2  = np.nan_to_num(st2.data.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
#    # Build resultant horizontal component
#    res_data = np.sqrt(data1 **2  + data2**2)
#    # Create a new Trace
#    tr = Trace(data=res_data, header=st1.stats)
#    
#    # Put it into a Stream
#    st = Stream([tr])
#    inv= read_inventory(xmlfile)
#    #Compute the ppsd
#    ppsd_z = PPSD(stZ[0].stats, metadata=inv, ppsd_length=3600)
#    ppsd_1 = PPSD(st[0].stats, metadata=inv, ppsd_length=3600)
#    #ppsd_1 = PPSD(st1[0].stats, metadata=inv, ppsd_length=3600)
#    #ppsd_2 = PPSD(st2[0].stats, metadata=inv, ppsd_length=3600)
#    #get the periods
#    # 3. Add data to PPSD
#    # ---------------------------------------------------------
#    ppsd_z.add(stZ)
#    ppsd_1.add(st1)
#    #ppsd_2.add(st2)
#    #get the mean of the periods and mean of the psd and plot them
#    T_z, mean_psd_z = ppsd_z.get_mean()
#    T_1, mean_psd_1 = ppsd_1.get_mean()
#    #T_2, mean_psd_2 = ppsd_2.get_mean()
#
#    # Convert from dB to linear
#    #mean_psd_z = 10 ** (mean_psd_z / 10)
#    #mean_psd_1 = 10 ** (mean_psd_1/ 10)
#    #mean_psd_2 = 10 ** (mean_psd_2 / 10)
#    # 4. Calculate H/V Ratio
#    # Standard formula: H = sqrt(N^2 + E^2) 
#    # Here we use the root mean square of horizontals
#    #h_mean    = np.sqrt(mean_psd_1**2 + mean_psd_2**2)
#    #hv_ratio  = h_mean / mean_psd_z
#    hv_ratio  = mean_psd_1 / mean_psd_z
#    # hv_ratio  = h_mean / np.abs(mean_psd_z)
#    #del st
#    #del res_data
#    #---------------------------------------------------------
#    ##return the result
#    return(T_z, hv_ratio)

def ppsd_style_smoothing(freq, amp, bins_per_octave=12):
    """
    Smooth amplitude spectrum using log-frequency binning
    similar to ObsPy PPSD behavior.
    """
    # Avoid zeros
    amp = np.maximum(amp, 1e-20)

    # Remove f=0
    freq = freq[freq > 0]
    amp = amp[:len(freq)]

    # Log-frequency axis
    logf = np.log10(freq)
    logf_min, logf_max = logf.min(), logf.max()

    # Number of bins
    n_bins = int((logf_max - logf_min) * bins_per_octave)

    # Log-spaced bin edges
    bin_edges = np.logspace(logf_min, logf_max, n_bins)

    freq_smooth = []
    amp_smooth = []

    # Bin and average
    for i in range(len(bin_edges) - 1):
        mask = (freq >= bin_edges[i]) & (freq < bin_edges[i+1])
        if np.any(mask):
            freq_smooth.append(np.mean(freq[mask]))
            amp_smooth.append(np.mean(amp[mask]))

    return np.array(freq_smooth), np.array(amp_smooth)

def smooth_log_freq(freq, amp, window=400):
    log_amp = np.log10(amp + 1e-20)
    smooth  = np.convolve(log_amp, np.ones(window)/window, mode='same')
    #smooth = savgol_filter(log_amp, window_length=31, polyorder=3)
    return 10**smooth




def calculate_hvr_psd(file_Z, file_1, file_2, xmlfile, amplification=None):
    #read the mseedfile
    stZ= read(file_Z)
    st1= read(file_1)
    st2= read(file_2)
    # Put it into a Stream
    #st = Stream([tr])
    inv= read_inventory(xmlfile)
    #Compute the ppsd
    ppsd_z = PPSD(stZ[0].stats, metadata=inv, ppsd_length=3600)
    ppsd_1 = PPSD(st1[0].stats, metadata=inv, ppsd_length=3600)
    ppsd_2 = PPSD(st2[0].stats, metadata=inv, ppsd_length=3600)
    #get the periods
    # 3. Add data to PPSD
    # ---------------------------------------------------------
    ppsd_z.add(stZ)
    ppsd_1.add(st1)
    ppsd_2.add(st2)
    #get the mean of the periods and mean of the psd and plot them
    T_z, mean_psd_z = ppsd_z.get_mean()
    T_1, mean_psd_1 = ppsd_1.get_mean()
    T_2, mean_psd_2 = ppsd_2.get_mean()

   # # Convert from dB to linear
    # Convert from dB to linear
    mean_psd_z = 10 ** (mean_psd_z / 10)
    mean_psd_1 = 10 ** (mean_psd_1 / 10)
    mean_psd_2 = 10 ** (mean_psd_2 / 10)

    # Clean NaNs
    #mean_psd_z = np.nan_to_num(mean_psd_z, nan=1e-20)
    #mean_psd_1 = np.nan_to_num(mean_psd_1, nan=0.0)
    #mean_psd_2 = np.nan_to_num(mean_psd_2, nan=0.0)
    # Compute H/V
    ## 4. Calculate H/V Ratio
    #h_mean = np.sqrt((mean_psd_1 + mean_psd_2) / 2)
    h_mean    = np.sqrt(mean_psd_1**2 + mean_psd_2**2)
    hv_ratio = h_mean / mean_psd_z
    #compute the amplification---------------------------------------------------------
    if(amplification):
        ##Resampling the data
        tr_a        =  st1[0].resample(sampling_rate =10.0, window = 'hann', no_filter = True)
        # --- Amplitude spectrum of Z ---
        n         = len(tr_a.data)
        dt        = tr_a.stats.delta
        freq_amp  = np.fft.rfftfreq(n, dt)
        ampZ      = np.abs(np.fft.rfft(tr_a.data))
        ################
        ##ampZ_norm = ampZ / ampZ.max()
        #ampZ_norm = konno_ohmachi_smoothing(ampZ, freq_amp, bandwidth=40)
        ampZ_norm = smooth_log_freq(freq_amp, ampZ)
        #freq_amp, ampZ_norm =    ppsd_style_smoothing(freq_amp, ampZ, bins_per_octave=12)
        mask      = freq_amp > 0
        T_amp     = 1./freq_amp[mask]
        ampZ_norm = ampZ_norm[mask]
        #print(T_amp)
        #print('-----' * 30)
        #print(ampZ_norm)
        #print(len(T_amp), len(ampZ_norm))
        #exit()
        #mask     = T < 4        # period < 4s
        #mask     = (T >= 1e-1) & (T < 4)        # period < 4s
        #T_sel    = T[mask]
        #HVR_sel  = HVR[mask]
        return (T_1, hv_ratio, T_amp, ampZ_norm)
    else:
        ##return the result
        return(T_1, hv_ratio)

def calculate_coherence(tr1, tr2):
    # ---------------------------------------------------------
    fs    = 250  # adjust depending on your sampling rate
    T_max = 1000  #second
    #Compute the nperseg as function sampling rate and the maximum periods
    nperseg = int(fs * T_max)
    # ---------------------------------------------------------
    # 3. Compute Horizontal Spectrum (root-mean-square)
    # Compute coherence
    f, Cxy = signal.coherence(tr1, tr2, fs, nperseg=nperseg)
    #f, Cxy = signal.coherence(trZ, trP, fs, nperseg=nperseg, nfft=3600)
    ######################################################
    # Remove f = 0 (avoid division by zero)
    mask = f > 0
    f    = f[mask]
    Cxy  = Cxy[mask]
    # Convert frequency to period
    T    = 1.0 / f
    ## Select period range between 1e-3 and 1000 seconds
    #mask_period = (T >= 1e-3) & (T <= 1000)
    #T           = T[mask_period]
    #Cxy         = Cxy[mask_period]
    return(T, Cxy)

########################################################################
#Set a component you need
component_Z      = "HH1" #Vertical compenent of the seismogram
compoent_2     = "HH2" #Vertical compenent of the seismogram
compoent_3     = "HH3" #Vertical compenent of the seismogram
#DEPTH of the seismic station in meter (m) in in the Lake Lucerne
#STATIONS_DEPTH = {"MUS03": 59, "MUS06": 71, "MUS07": 69, "MUS05": 49, "MUS04":41, "MUS08": 62, "MUS09": 59, "MUS10": 73}
STATIONS_DEPTH = {"MUS03": 59, "MUS06": 71, "MUS08": 62}
#STATIONS_DEPTH = {"MUS07": 69, "MUS05": 49, "MUS04":41, "MUS09": 59, "MUS10": 73}
#open a dictionary
#FILEDICT  = dict()
#path      = Pathdir("SeismicLucerne")
#
#for item in path.iterdir():
#    item_new = "%s/mseed"%(item)
#    #print(item_new)
#    path_years = Pathdir(item_new)
#    for item_year in path_years.iterdir():
#        #print(item_year)
#        path_months = Pathdir(item_year)
#        for item_month in path_months.iterdir():
#            #print(item_month)
#            path_days = Pathdir(item_month)
#            for item_day in path_days.iterdir():
#                #print(item_day)
#                for File_Z in Pathdir(item_day).glob("*%s*.mseed"%(component_Z)):
#                    basename = os.path.basename(str(File_Z))
#                    parts    = basename.split("_")
#                    station  = parts[1]
#                    File_2   = str(File_Z).replace('%s_%s'%(station,component_Z), '%s_%s'%(station,compoent_2)) 
#                    File_3   = str(File_Z).replace('%s_%s'%(station,component_Z), '%s_%s'%(station,compoent_3)) 
#                    #File_new = str(File).replace(component, compoent_P) 
#                    if(station in STATIONS_DEPTH):
#                        if station not in FILEDICT:
#                            FILEDICT[station] = []
#                        FILEDICT[station].append((str(File_Z), str(File_2), str(File_3)))
#




FILEDICT = {"MUS03": [('SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH1-2023-10-23.mseed',
                      'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed',
                      'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed')]}



#set params values
xlims  = (0.07, 10)        # works perfectly on log scale
ymin   = 0        # works perfectly on log scale
#plt.xlabel('frequency [Hz]')
xlabel = 'Period (s)'
#plt.ylabel('Coherence')
ylabel = 'HVSR (V/H)'
####################
nfigs  = 1
######################
fig_size= (10, 8)

# ---------------------------------------------------------
# 1. Load your three-component seismic data
#xmlfile
xmlfile = "RESPOSE-FILE.XML/xj-fdsn.xml"
# ---------------------------------------------------------
# Replace with your own file or stream

# Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"
#######################

for sta in FILEDICT:
    FILES_LIST_S   = [item[0] for item in FILEDICT[sta]] 
    #print(FILES_LIST_S)
    #exit()
    dates = [re.search(date_pattern, str(path)).group() for path in FILES_LIST_S if path and re.search(date_pattern, str(path))]
    ######################
    dates = sorted(set(dates))
    #####set the figure file name##################################
    figname = "Fig_HVR_%s_TEST.pdf"%(sta)
    #create a figure
    ##################Create a figure#######
    fig, ax  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
    ax_1     = ax.twinx()
    ######################################
    for item in FILEDICT[sta]:
        FILE_z = item[0]
        FILE_1 = item[1]
        FILE_2 = item[2]
        ################################
        try:
            #2.Run the calculation
            #freq, HVSR = calculate_hvr_psd(FILE_z, FILE_1, FILE_2, xmlfile)
            freq, HVSR, freq_a, amp_ = calculate_hvr_psd(FILE_z, FILE_1, FILE_2, xmlfile, amplification= True)
        except:
            continue
        # ---------------------------------------------------------
        #find the resonnance frequency
        f0, amp = find_resonance_frequency(freq, HVSR) 
        # ---- Quality control filters ----
        # --- QUALITY CONTROL BLOCK ---
        #Skip if amp is NaN, too high, or too low
#        if np.isnan(amp) or amp > 15.0 or amp < 1.1:
#            continue
#
#        #Optional: Skip if the curve is too "jagged"
#        if np.std(np.gradient(HVSR)) > 1.5:
#            continue
#        #detect spikes
#        spikes = np.sum(np.abs(np.diff(HVSR)) > 5)
#        if spikes > 10:
#            continue
#        #5. Plot HVSR curve
        #Plot the data on the figure
        #ax.semilogx(freq, HVSR, linewidth=1)
        ax.loglog(freq, HVSR, linewidth=1, color="darkblue")
        #Plot the resonnence frequency
        #print(len(freq_a), len(amp_))
        #exit()
        ax_1.loglog(freq_a, amp_, linewidth=1, color="r")
        #ax.axvline(x=f0, ymax=amp, color='grey',  linewidth=0.2, alpha=0.3)
        #ax.axvline(x=f0, ymax=amp,  linewidth=0.2, alpha=0.3)
        ax.axvline(x=f0, linewidth=0.4)
    ####################################################
    #ax.set_title('%s: %s--%s'%(sta, dates[0], dates[1]),fontsize=Figsize)
    ax.set_title('%s '%(sta),fontsize=Figsize)
    #Plot the figure
    #Make the grid of the axis
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    #Set ylim of x-axis
    ax.set_xlim(xlims[0], xlims[1])
#    ax.set_ylim(bottom = ymin)
    #disable ticks on  x-axis
    #Remove the ticks marks (small vertical lines) on the x-axis
    #ax.tick_params(axis="both", length = 0, color= "white", width = 0)
    ax.tick_params(axis='both', labelsize=Figsize)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    ax.set_xlabel(xlabel, labelpad = pad, fontsize = Figsize)
    ##########################################################
    #Save the figure
    plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
