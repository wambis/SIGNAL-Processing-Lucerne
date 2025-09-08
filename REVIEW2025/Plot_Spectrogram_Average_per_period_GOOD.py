#Plot the mseed seismic data files

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import MaxNLocator
from obspy import read, read_inventory
import matplotlib.pyplot as plt
from obspy import read, read_inventory
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy import read, read_inventory
from obspy.signal import PPSD
import io
from datetime import datetime, timedelta
from scipy.interpolate import interp1d





def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations

#Load 24-hour waveform data (MiniSEED)
#waveform = read("example.mseed")
waveform = read("DATASELECT/STATION03/XJ_MUS03_HH2-2023-10-23.mseed")

#Load station metadata (StationXML file)
inventory = read_inventory("RESPOSE-FILE.XML/xj-fdsn.xml")

#Select first trace
tr = waveform[0]

#Create PPSD object
ppsd = PPSD(tr.stats, inventory)

#Add entire 24-hour waveform data
ppsd.add(waveform)

#Define time extent for 24-hour spectrogram
# Define your 24-hour window

# Original data
#psd_matrix = np.array(ppsd.psd_values)       # shape: (n_times, n_periods)
#periods    = np.array(ppsd.period_bin_centers)
#times      = [t.datetime for t in ppsd.times_processed]
#time_nums  = mdates.date2num(times)


#SELECTED FIEXED PERIODs
target_periods = [100]




periods, mean_psd = ppsd.get_mean() 

#compute the average for periods > 10 seconds
# Create a mask for periods > 10 seconds
mask         = periods > 10.0
# Apply mask to periods and corresponding mean PSD
long_periods = periods[mask]
long_psd     = mean_psd[mask]
#Compute overall average PSD in that range
average_psd_over_10s = np.mean(long_psd)


#loop over the selected periods
for target_period in target_periods:
    #create a new matrix
    # Find index of closest period
    idx = (np.abs(periods - target_period)).argmin()
    #Get corresponding PSD value
    psd_at_target_period = mean_psd[idx]
    print(psd_at_target_period)
    #set the file name
    figname = "spectrogram_PSD_Meam_%ss.png"%(target_period) 
    #Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 5))
    #Set custom major ticks every 2 hours
    #plot the averge
    ax.semilogx(periods, mean_psd)
    #ax.axvline(periods[idx], color='red', linestyle='--', label=f"{periods[idx]:.2f} s")
    x1= 90; x2= 100
    ax.axvspan(x1, x2, facecolor="r", alpha = 0.2)
    #######################################
    ax.axhline(average_psd_over_10s, color='red', linestyle='--', label=f"Avg: {average_psd_over_10s:.2f} dB")
    #write the value on the y-axis
    # Add label at the left edge (x=0 in axis coordinates)
    ax.text(
        0, average_psd_over_10s +2,  # (x in axis coords, y in data coords)
        f"{average_psd_over_10s:.2f} dB",
        va='center', ha='left',
        fontsize=10, color='red',
        transform=ax.get_yaxis_transform()
    )
    #ax.xaxis.set_major_locator(locator)
    #ax.xaxis.set_major_formatter(formatter)
    ax.grid(True, which='both')
    #Label the Y-axis properly
    #set xlim
    ax.set_xlim(1, )
    #ax.set_xlim(1, 300 )
    #ax.set_yscale('log')
    ax.set_xlabel("Period (s)")
    #plt.ylabel("Mean PSD [dB]")
    plt.ylabel("PSD [dB ref m²/s⁴/Hz]")
    #####################################
    #Save the figure with high resolution
    fig.savefig(figname, dpi=300, bbox_inches="tight")
