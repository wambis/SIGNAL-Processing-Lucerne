##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
from scipy.signal import welch
from obspy import read
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
    # Select period range between 1e-3 and 1000 seconds
    mask_period = (T >= 1e-3) & (T <= 1000)
    T           = T[mask_period]
    Cxy         = Cxy[mask_period]
    return(T, Cxy)

########################################################################
#Set a component you need
component      = "HH1" #Vertical compenent of the seismogram
compoent_P     = "HDH" #Vertical compenent of the seismogram
#DEPTH of the seismic station in meter (m) in in the Lake Lucerne
#STATIONS_DEPTH = {"MUS03": 59, "MUS06": 71, "MUS07": 69, "MUS05": 49, "MUS04":41, "MUS08": 62, "MUS09": 59, "MUS10": 73}
STATIONS_DEPTH = {"MUS07": 69, "MUS05": 49, "MUS04":41, "MUS09": 59, "MUS10": 73}
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
                for File in Pathdir(item_day).glob("*%s*.mseed"%(component)):
                    basename = os.path.basename(str(File))
                    parts    = basename.split("_")
                    station  = parts[1]
                    File_new = str(File).replace('%s_%s'%(station,component), '%s_%s'%(station,compoent_P)) 
                    #File_new = str(File).replace(component, compoent_P) 
                    if(station in STATIONS_DEPTH):
                        if station not in FILEDICT:
                            FILEDICT[station] = []
                        FILEDICT[station].append((str(File), str(File_new)))



#print(FILEDICT)

#set params values
xlims  = (0.2, 200)        # works perfectly on log scale
#ymin   = 1e-3        # works perfectly on log scale
ymin   = 1e-4        # works perfectly on log scale
#plt.xlabel('frequency [Hz]')
xlabel = 'Period (s)'
#plt.ylabel('Coherence')
ylabel = 'P-Z Coherence'
#ylabel = 'X-Y Coherence'
#ylabel = 'P-X Coherence'
####################
nfigs  = 1
######################
fig_size= (10, 8)

# ---------------------------------------------------------
# 1. Load your three-component seismic data
# ---------------------------------------------------------
# Replace with your own file or stream

# Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"
#######################

for sta in FILEDICT:
    FILES_LIST_S   = [item[0] for item in FILEDICT[sta]] 
    #FILES_LIST_P   = [item[1] for item in FILEDICT[sta]] 
    #print(FILES_LIST_P)
    #exit()
    dates = [re.search(date_pattern, str(path)).group() for path in FILES_LIST_S if path and re.search(date_pattern, str(path))]
    ######################
    dates = sorted(set(dates))
    #####set the figure file name##################################
    figname = "Fig_coh_%s-%sx%s-%s_%s.pdf"%(sta, component,compoent_P, dates[0], dates[-1])
    #create a figure
    ##################Create a figure#######
    fig, ax  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
    ######################################
    for item in FILEDICT[sta]:
        FILE_S = item[0]
        FILE_P = item[1]

        try:
            tr_s = load_data(FILE_S)
            tr_p = load_data(FILE_P)
        except:
            continue
        #2.Run the calculation
        Periods, coh = calculate_coherence(tr_s, tr_p)
        #Plot the data on the figure
        #ax.loglog(T, Cxy, linewidth=0.1, color="darkblue")
        ax.loglog(Periods, coh, linewidth=0.1)
    ####################################################
    ax.set_title('%s: %s--%s'%(sta, dates[0], dates[1]),fontsize=Figsize)
    #Plot the figure
    #Make the grid of the axis
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    #Set ylim of x-axis
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(bottom = ymin)
    #disable ticks on  x-axis
    #Remove the ticks marks (small vertical lines) on the x-axis
    #ax.tick_params(axis="both", length = 0, color= "white", width = 0)
    ax.tick_params(axis='both', labelsize=Figsize)
    ax.set_ylabel(ylabel, labelpad = pad, fontsize = Figsize)
    ax.set_xlabel(xlabel, labelpad = pad, fontsize = Figsize)
    ##########################################################
    #Save the figure
    plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
