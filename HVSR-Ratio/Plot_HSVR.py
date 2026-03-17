##############################
#from hvseispy.hvsr import HVSR
import matplotlib.pyplot as plt
#from hvseispy.IO import read_sac, read_mseed, read_file
#from hvseispy.seismic import hvsr
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from obspy import read, Stream
#from hvseispy.plot_tools import plot_signal, plot_fft, plot_hv

import pathlib

import numpy as np
import matplotlib.pyplot as plt

import hvsrpy

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)




# Input file name (may be a relative or full path).
fnames = [["./XJ_MUS03_2023-10-23.mseed"]] 

print(f"Number of recordings: {len(fnames)}")
for fname_set in fnames:
    for file in fname_set:
        if not pathlib.Path(file).exists():
            raise FileNotFoundError(f"file {file} not found; check spelling.")
print("All files exist.")

#Preprocessing Settings 

preprocessing_settings = hvsrpy.settings.HvsrPreProcessingSettings()
preprocessing_settings.detrend = "linear"
preprocessing_settings.window_length_in_seconds = 100
preprocessing_settings.orient_to_degrees_from_north = 0.0
preprocessing_settings.filter_corner_frequencies_in_hz = (None, None)
preprocessing_settings.ignore_dissimilar_time_step_warning = False

print("Preprocessing Summary")
print("-"*60)
preprocessing_settings.psummary() 

#Processing Settings

processing_settings = hvsrpy.settings.HvsrTraditionalProcessingSettings()
processing_settings.window_type_and_width = ("tukey", 0.2)
processing_settings.smoothing=dict(operator="konno_and_ohmachi",
                                   bandwidth=40,
                                   center_frequencies_in_hz=np.geomspace(0.2, 50, 200))
processing_settings.method_to_combine_horizontals = "geometric_mean"
processing_settings.handle_dissimilar_time_steps_by = "frequency_domain_resampling"

print("Processing Summary")
print("-"*60)
processing_settings.psummary() 

#Calculate HVSR


srecords = hvsrpy.read(fnames)
srecords = hvsrpy.preprocess(srecords, preprocessing_settings)
hvsr = hvsrpy.process(srecords, processing_settings)


print("\nStatistical Summary:")
print("-"*20)
hvsrpy.summarize_hvsr_statistics(hvsr)
(fig, ax) = hvsrpy.plot_single_panel_hvsr_curves(hvsr,)
ax.get_legend().remove()
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

# Save the figure
fig.savefig("hvsr_result.png", dpi=300, bbox_inches="tight")

#plt.show() 
print("wamba")
exit()

# Example: Reading MSEED files

## Example: Reading MSEED files
#north_    = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH2-2023-10-23.mseed'
#east_     = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH3-2023-10-23.mseed'
#vertical_ = 'SeismicLucerne/muota-SingSta-03/mseed/2023/10/23/XJ_MUS03_HH1-2023-10-23.mseed'
#
#
#
#trZ = read(vertical_)[0]
#trN = read(north_)[0]
#trE = read(east_)[0]
#
##channels = [tr.stats.channel for tr in st]
#
#trZ.stats.channel = "Z"
#trN.stats.channel = "N"
#trE.stats.channel = "E"
#
#st = Stream([trN, trE, trZ]) 
#
#
##write into a file
#st.write("XJ_MUS03_2023-10-23.mseed", format="MSEED")
#
#
#
## Example: Reading MSEED files
#north, vertical, east = read_mseed('XJ_MUS03_2023-10-23.mseed')

file_mseed = "XJ_MUS03_2023-10-23.mseed"

#print(trE.stats.sampling_rate)
#
#print(north.shape, vertical.shape, east.shape)
#print(len(east), len(vertical), len(north))
#print((len(north) - 1) / (north_time[-1] - north_time[0]))
#exit()

#dt = 0.004
dt = 1/324.0
win_len = 20.0
win_samples = int(win_len / dt)

# Trim to a multiple of window size
new_len = (len(north) // win_samples) * win_samples

north = np.asarray(north[:new_len], dtype=float)
vertical = np.asarray(vertical[:new_len], dtype=float)
east = np.asarray(east[:new_len], dtype=float)

hv_mean, hv, freq = hvsr(
    acc_data=[north, vertical, east],
    dt=dt,
    win_len=win_len,
    taper=('cosine', 0.05),
    smooth_bandwidth=40.0,
    overlap=0.5,
    fftmin=0.1,
    fftmax=50.0,
    hvmin=0.1,
    hvmax=10.0
)





##print(north)
##print(vertical)
#print(type(east))
#print(north.shape, vertical.shape, east.shape)
##print(len(east), len(vertical), len(north))
#exit()
#
#
## Example: Calculating HVSR
#hv_mean, hv, freq = hvsr(
#    acc_data=[north, vertical, east],  # Input data as a list of numpy arrays
#    dt=0.004,  # Sampling period in seconds
#    win_len=20.0,  # Window length in seconds
#    taper=('cosine', 0.05),  # Taper type and amount
#    smooth_bandwidth=40.0,  # Smoothing bandwidth for Konno-Ohmachi smoothing
#    overlap=0.5,  # Overlap between windows (0 to 1)
#    fftmin=0.1,  # Minimum frequency for FFT
#    fftmax=50.0,  # Maximum frequency for FFT
#    hvmin=0.1,  # Minimum frequency for HVSR
#    hvmax=10.0  # Maximum frequency for HVSR
#) 
#
#
##print(len(east), len(vertical), len(north))
#print(hv_mean)
#exit()

fig = plot_hv(
    hv_mean, hv, freq,
    fmin=0.1, fmax=10.0,
    name='Station_Name',
    plot_windows=True
)

fig.savefig("hvsr_station_name.png", dpi=300, bbox_inches="tight")




## Step 1: Load seismic data
## Replace 'your_data.mseed' with your actual file path
#hvsr = HVSR()
#hvsr.load_data('your_data.mseed')
#
## Step 2: Preprocess and compute H/V ratio
#hvsr.preprocess_data()
#hvsr.compute_hvsr()
#
## Step 3: Plot the H/V spectral ratio
#fig, ax = plt.subplots(figsize=(10, 6))
#ax.plot(hvsr.freq, hvsr.hvsr, label='H/V Ratio', color='blue')
#ax.set_xscale('log')
#ax.set_xlabel('Frequency (Hz)')
#ax.set_ylabel('H/V Spectral Ratio')
#ax.set_title('Horizontal-to-Vertical Spectral Ratio')
#ax.grid(True, which='both', linestyle='--', linewidth=0.5)
#ax.legend()
#plt.show()
#
