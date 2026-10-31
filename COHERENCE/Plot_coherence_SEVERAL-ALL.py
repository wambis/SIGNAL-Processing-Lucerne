import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from scipy import signal
import os, re
from pathlib import Path as Pathdir
from datetime import datetime

#################################
Figsize = 14
fticks  = 18
pad     = 10

#Pattern for YYYY-MM-DD
date_pattern = r"\d{4}-\d{2}-\d{2}"

##########defined a function #################
def compute_coherence(tr1, tr2, fs):
    T_max = 1000
    nperseg = int(fs * T_max)

    f, Cxy = signal.coherence(tr1, tr2, fs, nperseg=nperseg)

    mask = f > 0
    f = f[mask]
    Cxy = Cxy[mask]

    T = 1.0 / f
    mask_period = (T >= 1e-3) & (T <= 1000)

    return T[mask_period], Cxy[mask_period]


def QUALITY_CONTROL(T, Cxy):
    mask = (T >= 8.5) & (T <= 50)
    return T[mask], Cxy[mask]


# -----------------------------
# Prepare figure
# -----------------------------
fig, ax = plt.subplots(figsize=(10, 8))

xlabel = "Period (s)"
#ylabel = "X–Z Coherence"
ylabel = "Z–P Coherence"

component_Z = "HH1"
component_2 = "HH2"
component_P = "HDH"
#set the priods band
Tmin = 5  #in second
Tmax = 500  #in second
#set the limit of the axis
xlims  = (5 , 500)
#ymin   = 1e-4
ymin   = 1e-2 * 2

ylims =(1e-2 * 2,  1)

path = Pathdir("SeismicLucerne")

# -----------------------------
# Storage for averaging
# -----------------------------
all_T = []
all_Cxy = []

dates = []
# -----------------------------
# Loop through all files
# -----------------------------
for item in path.iterdir():
    path_years = Pathdir(f"{item}/mseed")
    for item_year in path_years.iterdir():
        for item_month in Pathdir(item_year).iterdir():
            for item_day in Pathdir(item_month).iterdir():

                for File_Z in Pathdir(item_day).glob(f"*{component_Z}*.mseed"):

                    basename = os.path.basename(str(File_Z))
                    station  = basename.split("_")[1]

                    File_2 = str(File_Z).replace(f"{station}_{component_Z}",
                                                 f"{station}_{component_2}")
                    File_P = str(File_Z).replace(f"{station}_{component_Z}",
                                                 f"{station}_{component_P}")

                    trZ = read(File_Z)[0]
                    trN = read(File_2)[0]
                    trP = read(File_P)[0]

                    #filter the traces
                    trZ.detrend("linear"); trN.detrend("linear"); trP.detrend("linear")
                    trZ.taper(max_percentage=0.05); trN.taper(max_percentage=0.05); trP.taper(max_percentage=0.05)
                    ####Filter the traces #########################
                    trZ.filter("bandpass", freqmin=1./Tmax, freqmax=1/Tmin)
                    trN.filter("bandpass", freqmin=1./Tmax, freqmax=1/Tmin)
                    trP.filter("bandpass", freqmin=1./Tmax, freqmax=1/Tmin)

                    fs = trZ.stats.sampling_rate

                    ## Compute coherence X-Z
                    #T, Cxy = compute_coherence(trN, trZ, fs)

                    # Compute coherence X-Z
                    T, Cxy = compute_coherence(trZ, trP, fs)

                    #Store for averaging
                    all_T.append(T)
                    all_Cxy.append(Cxy)

                    #grab the date
                    basename = os.path.basename(File_Z)
                    #station  = basename.split("_")[1]
                    #######################
                    try:
                        date = [re.search(date_pattern, str(path)).group() for path in [basename] if path and re.search(date_pattern, str(path))]
                        dates.append(date[0])
                    except:
                        pass
                    #print(dates)
                    #exit()


# -----------------------------
# Interpolate all curves to a common T grid
# -----------------------------
T_common = np.logspace(np.log10(0.001), np.log10(1000), 2000)

Cxy_interp = []
for T, Cxy in zip(all_T, all_Cxy):
    Cxy_interp.append(np.interp(T_common, T[::-1], Cxy[::-1]))

Cxy_mean = np.mean(Cxy_interp, axis=0)

# -----------------------------
# Plot ONE final coherence curve
# -----------------------------
#ax.loglog(T_common, Cxy_mean, color="k", linewidth=2.0, label="Mean coherence")
ax.loglog(T_common, Cxy_mean, color="k", linewidth=2.0)

## Highlight high-coherence region
#T_sel, Cxy_sel = QUALITY_CONTROL(T_common, Cxy_mean)
#ax.loglog(T_sel, Cxy_sel, color="r", linewidth=2.5, label="High coherence")

# -----------------------------
# Final formatting
# -----------------------------
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.set_xlim(xlims)
#ax.set_ylim(bottom=ymin)
ax.set_ylim(ylims)

#order the dates
# Convert to datetime objects and sort
dates = sorted(dates, key=lambda d: datetime.strptime(d, "%Y-%m-%d"))
#set title
ax.set_title("%s—%s"%(dates[0], dates[-1]), fontsize=Figsize)

ax.set_xlabel(xlabel, fontsize=fticks, labelpad=pad)
ax.set_ylabel(ylabel, fontsize=fticks, labelpad=pad)
ax.tick_params(axis="both", labelsize=fticks)

#ax.legend(fontsize=fticks)

plt.tight_layout()
plt.savefig("Coherence-AVERAGED.pdf", dpi=300, bbox_inches="tight")
#plt.show()
