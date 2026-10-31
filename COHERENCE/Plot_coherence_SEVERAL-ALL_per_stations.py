import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from scipy import signal
import os, re
from pathlib import Path as Pathdir
from datetime import datetime
#######################################
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

date_pattern = r"\d{4}-\d{2}-\d{2}"

#################################
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

#################################
# PARAMETERS
#################################
xlabel = "Period (s)"
ylabel = "Z–P Coherence"
#ylabel = "X–Z Coherence"

component_Z = "HH1"
component_2 = "HH2"
component_P = "HDH"

Tmin = 5
Tmax = 500
xlims = (5, 500)
ylims = (1e-3 * 6, 1)

path = Pathdir("SeismicLucerne")

#################################
# STORAGE PER STATION
#################################
#stations = ["MUS03", "MUS04", "MUS05", "MUS06", "MUS07","MUS08", "MUS09", "MUS10"]
stations = ["MUS03", "MUS04", "MUS05", "MUS06","MUS08", "MUS10"]

station_T = {sta: [] for sta in stations}
station_Cxy = {sta: [] for sta in stations}

dates = []

#################################
# LOOP THROUGH FILES
#################################
for item in path.iterdir():
    path_years = Pathdir(f"{item}/mseed")
    for item_year in path_years.iterdir():
        for item_month in Pathdir(item_year).iterdir():
            for item_day in Pathdir(item_month).iterdir():

                for File_Z in Pathdir(item_day).glob(f"*{component_Z}*.mseed"):

                    basename = os.path.basename(str(File_Z))
                    station  = basename.split("_")[1]
                    #print(station)
                    #exit()

                    if station not in stations:
                        continue

                    File_2 = str(File_Z).replace(f"{station}_{component_Z}",
                                                 f"{station}_{component_2}")
                    File_P = str(File_Z).replace(f"{station}_{component_Z}",
                                                 f"{station}_{component_P}")

                    trZ = read(File_Z)[0]
                    trN = read(File_2)[0]
                    trP = read(File_P)[0]

                    #Filter
                    for tr in [trZ, trN, trP]:
                        tr.detrend("linear")
                        tr.taper(max_percentage=0.05)
                        tr.filter("bandpass", freqmin=1./Tmax, freqmax=1./Tmin)

                    fs = trZ.stats.sampling_rate

                    # Compute coherence Z–P
                    T, Cxy = compute_coherence(trZ, trP, fs)
                    #T, Cxy = compute_coherence(trZ, trN, fs)

                    station_T[station].append(T)
                    station_Cxy[station].append(Cxy)

                    # Extract date
                    try:
                        date = re.search(date_pattern, basename).group()
                        dates.append(date)
                    except:
                        pass

#################################
# SORT DATES
#################################
dates = sorted(dates, key=lambda d: datetime.strptime(d, "%Y-%m-%d"))

#################################
# COMMON PERIOD GRID
#################################
T_common = np.logspace(np.log10(5), np.log10(500), 1500)

#################################
# PLOT PER STATION
#################################
fig, ax = plt.subplots(figsize=(12, 8))

#colors = ["red", "blue","k", "c", "m", "green", "orange", "purple"]
colors = ["red", "blue","k", "c", "green", "orange"]


for i, sta in enumerate(stations):

    curves_T = station_T[sta]
    curves_C = station_Cxy[sta]

    if len(curves_T) == 0:
        continue

    Cxy_interp = []

    for T, Cxy in zip(curves_T, curves_C):

        # Sort T ascending
        idx = np.argsort(T)
        T_sorted = T[idx]
        Cxy_sorted = Cxy[idx]

        # Overlap mask
        mask = (T_sorted >= T_common.min()) & (T_sorted <= T_common.max())
        if np.sum(mask) < 10:
            continue

        T_use = T_sorted[mask]
        Cxy_use = Cxy_sorted[mask]

        Cxy_i = np.interp(T_common, T_use, Cxy_use)
        Cxy_interp.append(Cxy_i)

    if len(Cxy_interp) == 0:
        continue

    Cxy_mean = np.mean(Cxy_interp, axis=0)

    ax.loglog(T_common, Cxy_mean, color=colors[i], linewidth=1.2, label=sta)
    #ax.loglog(T_common, Cxy_mean, linewidth=1.0, label=sta)

#################################
# FINAL FORMATTING
#################################
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.set_xlim(xlims)
ax.set_ylim(ylims)

ax.set_title(f"{dates[0]} — {dates[-1]}", fontsize=Figsize)
ax.set_xlabel(xlabel, fontsize=fticks, labelpad=pad)
ax.set_ylabel(ylabel, fontsize=fticks, labelpad=pad)
ax.tick_params(axis="both", labelsize=fticks)

ax.legend(fontsize=fticks)

plt.tight_layout()
plt.savefig("Coherence-PER-STATION.pdf", dpi=300, bbox_inches="tight")
