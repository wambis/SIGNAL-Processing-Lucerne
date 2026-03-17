##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from obspy import read



def find_resonance_frequency(f, hvsr):

    idx = np.argmax(hvsr)
    f0 = f[idx]
    amp = hvsr[idx]

    return f0, amp

def calculate_hvr(trZ, tr1, tr2, periode=None):
    #Sampling frequency
    fs = trZ.stats.sampling_rate
    #T_max = 1000
    # ---------------------------------------------------------
    # 2. Compute PSD using Welch method (Nakamura requires PSD)
    # ---------------------------------------------------------
    #nperseg = 4096  # adjust depending on your sampling rate
    #nperseg = 2048  # adjust depending on your sampling rate
    #nperseg = 3600  # adjust depending on your sampling rate
    nperseg = 1024  # adjust depending on your sampling rate
    #nperseg = int(fs * T_max)  # adjust depending on your sampling rate
    f, Pz = welch(trZ.data, fs=fs, nperseg=nperseg, average='median')
    _, Pn = welch(tr1.data, fs=fs, nperseg=nperseg, average='median')
    _, Pe = welch(tr2.data, fs=fs, nperseg=nperseg, average='median')
    
    #---------------------------------------------------------
    # 3. Compute Horizontal Spectrum (root-mean-square)
    #---------------------------------------------------------
    Ph = np.sqrt((Pn + Pe)/2)
    #---------------------------------------------------------
    # 4. Compute HVSR (Nakamura)
    #---------------------------------------------------------
    HVSR = Ph / np.sqrt(Pz)
    if(periode):
        # Remove f = 0 (avoid division by zero)
        mask = f > 0
        f    = f[mask]
        HVSR = HVSR[mask]
        #Convert frequency to period
        T = 1.0 / f
        # Select period range between 1e-3 and 1000 seconds
        mask_period = (T >= 1e-3) & (T <= 1000)
        T    = T[mask_period]
        HVSR = HVSR[mask_period]
        return(T, HVSR)
    ##return the result
    return(f, HVSR)




#Example: Reading MSEED files
vertical = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH1-2023-10-23.mseed'
north    = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed'
east     = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed'

# ---------------------------------------------------------
# 1. Load your three-component seismic data
# ---------------------------------------------------------
# Replace with your own file or stream
#st = read("your_3C_file.mseed")

# Extract components (Z, N, E)
trZ = read(vertical)[0]
trN = read(north)[0]
trE = read(east)[0]






freq, HVSR = calculate_hvr(trZ, trN, trE, periode=True)
# ---------------------------------------------------------
#find the resonnance frequency
f0, amp = find_resonance_frequency(freq, HVSR) 
# 5. Plot HVSR curve
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))
#plt.semilogx(freq, HVSR, linewidth=2, color="darkblue")
plt.loglog(freq, HVSR, linewidth=2, color="darkblue")

#Plot the resonnence frequency
plt.axvline(x=f0, ymax=amp, color='grey',  linewidth=4, alpha=0.3)

#plt.xlabel("Frequency (Hz)")
plt.xlabel("Period (s)")
plt.ylabel("HVSR (H/V)")
plt.title("Horizontal-to-Vertical Spectral Ratio (Nakamura Method)")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
figname = "HVSR_GOOD.png"
#Save the figure
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
#plt.show()
