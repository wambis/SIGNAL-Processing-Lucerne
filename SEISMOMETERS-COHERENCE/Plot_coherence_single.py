##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
from scipy.signal import welch
from obspy import read
from scipy import signal
import os, re



#Set parameters
xy_coords   = (0.85, 0.955)
#xy_coords   = (0.85, 0.97)
#xy_coords_D = (0.55, 0.85)
#set the fontsize
Figsize= 18
pad    = 10
#define text size
f_size = 14
####for ticks #####
fticks = 24
#fticks = 18


def compute_coherence(tr1, tr2, fs):
    # ---------------------------------------------------------
    # 2. Compute PSD using Welch method (Nakamura requires PSD)
    # ---------------------------------------------------------
    #nperseg = 256  # adjust depending on your sampling rate
    #nperseg = 250000  # adjust depending on your sampling rate
    T_max = 1000  #second
    #Compute the nperseg as function sampling rate and the maximum periods
    nperseg = int(fs * T_max)
    
    # ---------------------------------------------------------
    # 3. Compute Horizontal Spectrum (root-mean-square)
    f, Cxy = signal.coherence(tr1, tr2, fs, nperseg=nperseg)
    #f, Cxy = signal.coherence(trZ, trP, fs, nperseg=nperseg, nfft=3600)
    ######################################################
    # Remove f = 0 (avoid division by zero)
    mask = f > 0
    f = f[mask]
    Cxy = Cxy[mask]
    
    # Convert frequency to period
    T = 1.0 / f
    
    # Select period range between 1e-3 and 1000 seconds
    mask_period = (T >= 1e-3) & (T <= 1000)
    T = T[mask_period]
    Cxy = Cxy[mask_period]
    return (T, Cxy)

#---------------------------------------------------------

def QUALITY_CONTROL(T, HVR):
    #Select only periods below 4 s
    #mask     = T < 4        # period < 4s
    mask     = (T >=10) & (T <=75)        # period < 4s
    #mask     = (T >=8.5) & (T <=50)        # period < 4s
    T_sel    = T[mask]
    HVR_sel  = HVR[mask]
    return(T_sel, HVR_sel)


#Example: Reading MSEED files
#north    = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed'
#east     = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed'

##Wind event
#FileH1 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH1-2024-03-29.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH2-2024-03-29.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HH3-2024-03-29.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-05/mseed/2024/03/29/XJ_MUS05_HDH-2024-03-29.mseed'


#FileH1 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH1-2024-01-13.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH2-2024-01-13.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HH3-2024-01-13.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-05/mseed/2024/01/13/XJ_MUS05_HDH-2024-01-13.mseed'

#FileH1 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH1-2023-09-22.mseed'
#FileH2 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH2-2023-09-22.mseed'
#FileH3 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HH3-2023-09-22.mseed'
#FileP  = 'SeismicLucerne/muota-SingSta-03/mseed/2023/09/22/XJ_MUS03_HDH-2023-09-22.mseed'

#Boats event
FileH1 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/02/XJ_MUS03_HH1-2023-10-02.mseed'
FileH2 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/02/XJ_MUS03_HH2-2023-10-02.mseed'
FileH3 = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/02/XJ_MUS03_HH3-2023-10-02.mseed'
FileP  = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/02/XJ_MUS03_HDH-2023-10-02.mseed'

# ---------------------------------------------------------
# 1. Load your three-component seismic data
# ---------------------------------------------------------
# Replace with your own file or stream

# Extract components (Z, N, E)
trZ = read(FileH1)[0]
trN = read(FileH2)[0]
trE = read(FileH3)[0]
trP = read(FileP)[0]

# Sampling frequency
fs = trZ.stats.sampling_rate


#Compute the cohernce
T, Cxy = compute_coherence(trZ, trP, fs)
#T, Cxy = compute_coherence(trN, trP, fs)
#T, Cxy = compute_coherence(trE, trZ, fs)
#T, Cxy = compute_coherence(trE, trN, fs)

#Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"
basename = os.path.basename(FileH1)
station  = basename.split("_")[1]

#######################
dates = [re.search(date_pattern, str(path)).group() for path in [basename] if path and re.search(date_pattern, str(path))]



#plot the high coherence
plot_high_coherence = True
#set params values
#xlims  = (0.1, 200)        # works perfectly on log scale
xlims  = (0.2, 200)        # works perfectly on log scale
#xlims  = (1, 200)        # works perfectly on log scale
#ymin   = 1e-3        # works perfectly on log scale
ymin   = 1e-4        # works perfectly on log scale
#plt.xlabel('frequency [Hz]')
xlabel = 'Period (s)'
#plt.ylabel('Coherence')
ylabel = 'P-Z Coherence'
#ylabel = 'X-Z Coherence'
#ylabel = 'P-X Coherence'
#ylabel = 'X-Y Coherence'
####################
nfigs  = 1
######################
fig_size= (10, 8)
##################Create a figure#######
fig, ax  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)


#set ylimit
#write the text on the figure
ax.annotate(
    station,     # Text for the annotation
    xy=  xy_coords,                   # Coordinates (relative to axes) for the annotation
    xycoords='axes fraction',    # Use axes fraction for positioning
    fontsize= Figsize                 # Font size for the annotation
    #bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
)
##Remove the frame on the plot
#ax.spines["top"].set_visible(False)
#ax.spines["right"].set_visible(False)
##################################
#Set_title
ax.set_title(dates[0],fontsize=Figsize)
#Plot the figure
#ax.loglog(T, Cxy, linewidth=1.0, color="darkblue")
ax.loglog(T, Cxy, linewidth=0.8, color="k")
#plot higher coherence
if(plot_high_coherence):
    T_sel, Cxy_sel = QUALITY_CONTROL(T, Cxy)
    #ax.loglog(T_sel, Cxy_sel, linewidth=1.0, color="r")
    #ax.loglog(T_sel, Cxy_sel, linewidth=1.0, color="k", alpha=0.5)
#number of vertical lines for grid
#Make the grid of the axis
#ax.grid(visible = True, axis = "x", alpha = 0.7)
ax.grid(True, which="both", linestyle="--", alpha=0.5)
#Set ylim of x-axis
ax.set_xlim(xlims[0], xlims[1])
ax.set_ylim(bottom = ymin)
#Set label
#disable ticks on  x-axis
#Remove the ticks marks (small vertical lines) on the x-axis
#ax.tick_params(axis="both", length = 0, color= "white", width = 0)
#ax.tick_params(axis='both', labelsize=Figsize)
ax.set_ylabel(ylabel, labelpad = pad, fontsize = fticks)
ax.set_xlabel(xlabel, labelpad = pad, fontsize = fticks)
##remove x-label
#ax.set_xticks([])
#ax.set_xticklabels([])
#Customize y-axis tick label font size
ax.tick_params(axis='both', labelsize=fticks)  # Change 14 to your desired size
#save the figure
figname = "Coherence-%s-%s-%s.pdf"%(station, dates[0], re.sub(r"\s+", "", ylabel))
#Save the figure
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
