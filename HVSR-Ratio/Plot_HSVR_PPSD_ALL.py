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
##############################
from concurrent.futures import ProcessPoolExecutor
import numpy as np

########################################################
AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}
#size for anotation
#xy_coords = (0.71, 0.95)
xy_coords   = (0.75, 0.95)
xy_coords_D = (0.75, 0.86)
#set the fontsize
Figsize= 18
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
#    st1= read(file_1)
#    st2= read(file_2)
#    inv= read_inventory(xmlfile)
#    #Compute the ppsd
#    ppsd_z = PPSD(stZ[0].stats, metadata=inv, ppsd_length=3600)
#    ppsd_1 = PPSD(st1[0].stats, metadata=inv, ppsd_length=3600)
#    ppsd_2 = PPSD(st2[0].stats, metadata=inv, ppsd_length=3600)
#    #get the periods
#    # 3. Add data to PPSD
#    # ---------------------------------------------------------
#    ppsd_z.add(stZ)
#    ppsd_1.add(st1)
#    ppsd_2.add(st2)
#    #get the mean of the periods and mean of the psd and plot them
#    T_z, mean_psd_z = ppsd_z.get_mean()
#    T_1, mean_psd_1 = ppsd_1.get_mean()
#    T_2, mean_psd_2 = ppsd_2.get_mean()
#
#    # Convert from dB to linear
#    mean_psd_z = 10 ** (mean_psd_z / 10)
#    mean_psd_1 = 10 ** (mean_psd_1/ 10)
#    mean_psd_2 = 10 ** (mean_psd_2 / 10)
#    # 4. Calculate H/V Ratio
#    # Standard formula: H = sqrt(N^2 + E^2) 
#    # Here we use the root mean square of horizontals
#    h_mean    = np.sqrt(mean_psd_1**2 + mean_psd_2**2)
#    hv_ratio  = h_mean / mean_psd_z
#    #---------------------------------------------------------
#    ##return the result
#    return(T_1, hv_ratio)


def calculate_hvr_psd_log(file_Z, file_1, file_2, xmlfile, periode=None):
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
    #---------------------------------------------------------
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


def QUALITY_CONTROL(T, HVR):
    #Select only periods below 4 s
    #mask     = T < 4        # period < 4s
    mask     = (T >= 1e-1) & (T < 4)        # period < 4s
    T_sel    = T[mask]
    HVR_sel  = HVR[mask]
    if(np.max(HVR_sel) >10.1):
        control = 1
    else:
        control = 0
    return(control, T_sel, HVR_sel)

########################################################################
#Set a component you need
component_Z      = "HH1" #Vertical compenent of the seismogram
compoent_2     = "HH2" #Vertical compenent of the seismogram
compoent_3     = "HH3" #Vertical compenent of the seismogram
#DEPTH of the seismic station in meter (m) in in the Lake Lucerne
STATIONS_DEPTH = {"MUS03": 59, "MUS06": 71, "MUS07": 69, "MUS05": 49, "MUS04":41, "MUS08": 62, "MUS09": 59, "MUS10": 73}
#STATIONS_DEPTH = {"MUS03": 59, "MUS05": 49, "MUS06": 71, "MUS08": 62}
#STATIONS_DEPTH = {"MUS07": 69, "MUS05": 49, "MUS04":41, "MUS09": 59, "MUS10": 73}
#open a dictionary
FILEDICT  = dict()
path      = Pathdir("SeismicLucerne")

for item in path.iterdir():
    item_new = "%s/mseed"%(item)
    #print(item_new)
    path_years = Pathdir(item_new)
    for item_year in path_years.iterdir():
        #print(item_year)
        path_months = Pathdir(item_year)
        for item_month in path_months.iterdir():
            #print(item_month)
            path_days = Pathdir(item_month)
            for item_day in path_days.iterdir():
                #print(item_day)
                for File_Z in Pathdir(item_day).glob("*%s*.mseed"%(component_Z)):
                    basename = os.path.basename(str(File_Z))
                    parts    = basename.split("_")
                    station  = parts[1]
                    File_2   = str(File_Z).replace('%s_%s'%(station,component_Z), '%s_%s'%(station,compoent_2)) 
                    File_3   = str(File_Z).replace('%s_%s'%(station,component_Z), '%s_%s'%(station,compoent_3)) 
                    #File_new = str(File).replace(component, compoent_P) 
                    if(station in STATIONS_DEPTH):
                        if station not in FILEDICT:
                            FILEDICT[station] = []
                        FILEDICT[station].append((str(File_Z), str(File_2), str(File_3)))




#set params values
#xlims  = (0.07, 10)        # works perfectly on log scale
xlims  = (0.07, 7)        # works perfectly on log scale
ymin   = 0        # works perfectly on log scale
ymax   = 100        # works perfectly on log scale
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
    figname = "Fig_HVR_%s_%s--%s.pdf"%(sta, dates[0], dates[-1])
    #create a figure
    ##################Create a figure#######
    fig, ax  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
    ######################################
    for item in FILEDICT[sta]:
        FILE_z = item[0]
        FILE_1 = item[1]
        FILE_2 = item[2]
        ################################
        try:
            #2.Run the calculation
            #freq, HVSR = calculate_hvr_psd(FILE_z, FILE_1, FILE_2, xmlfile)
            freq, HVSR = calculate_hvr_psd_log(FILE_z, FILE_1, FILE_2, xmlfile)
        except:
            continue
        # ---------------------------------------------------------
        # ---- Quality control filters ----
        # --- QUALITY CONTROL BLOCK ---
        c, freq_sel, HVR_sel = QUALITY_CONTROL(freq, HVSR)
        #Skip if amp is NaN, too high, or too low
        if(c==1):
            #find the resonnance frequency
            #f0, amp = find_resonance_frequency(freq_sel, HVR_sel) 
            # ---------------------------------------------------------
            #Plot the data on the figure
            #ax.semilogx(freq, HVSR, linewidth=1)
            ax.loglog(freq, HVSR, linestyle="--", linewidth=0.5)
            #ax.loglog(freq, HVSR, linewidth=1, color="darkblue")
            #Plot the resonnence frequency
            #ax.axvline(x=f0, ymax=amp, color='grey',  linewidth=0.2, alpha=0.3)
            #ax.axvline(x=f0, ymax=amp,  linewidth=0.2, alpha=0.3)
            #ax.axvline(x=f0, linewidth=0.4)
    ####################################################
    #ax.set_title('%s: %s--%s'%(sta, dates[0], dates[-1]),fontsize=Figsize)
    ax.set_title(' %s--%s'%(dates[0], dates[-1]),fontsize=Figsize)
    #Annotate the seismic station on the figure
    ax.annotate(
            sta,                        #Text for the annotation
            xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
            xycoords='axes fraction',         # Use axes fraction for positioning
            fontsize= f_size,                 # Font size for the annotation
            color   = 'k'  ,                  #text color
            alpha   = 0.9 )                   #alpha for the text
    #Plot the figure
    #Make the grid of the axis
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    #Set ylim of x-axis
    ax.set_xlim(xlims[0], xlims[1])
#    ax.set_ylim(bottom = ymin)
    ax.set_ylim(top = ymax)
    #disable ticks on  x-axis
    #Remove the ticks marks (small vertical lines) on the x-axis
    #ax.tick_params(axis="both", length = 0, color= "white", width = 0)
    ax.tick_params(axis='both', labelsize=Figsize)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    ax.set_xlabel(xlabel, labelpad = pad, fontsize = Figsize)
    ##########################################################
    #Save the figure
    plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
