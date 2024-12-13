import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdate
import re, pandas as pd 
import dolfyn as dlfn
import datetime

#
def verical_grid():
    #number of vertical lines
    num_vertical = 8
    #generate spaced x-values for vertical grid lines
    vertical_locations = plt.MaxNLocator(num_vertical +2)

    return vertical_locations


def start_endtime(time_serie):
    #minimum time
    tstart     = time_serie.min()
    #tstart     = time_serie[0]
    #Add one day to the minimum time
    tend      = tstart + pd.to_timedelta(1, unit= 'D')
    #tend       = time_serie[-1]
    #Convert to minimum
    tend       = tend.to_numpy()
    return (tstart, tend)



ADCP_FILE = "DATASELECT/MUO2_DOWN_11_23.000"


#Set the SEAFLOOR DEPTH at the location of the ADCP
###For the deploiement from `2023-09-16` to `2023-11-21` #########
SEAFLOOR_DEPTH= 120.0 # in meter
###Set the ADCP DEPTH
ADCP_UP_HEIGTH_FROM_LAKEBED= 20.9    #in meter
ADCP_DOWN_HEIGTH_FROM_LAKEBED= 20.0  #in meter


df    = dlfn.read(ADCP_FILE) 

target_day = "2023-10-23"


#get the backscatter
BS            = df.amp.sel(time = target_day)


#get the paramter
beam_angle   = df.beam_angle

blank_dist   = df.blank_dist

ranges       = BS.range.data
time         = BS.time.data
beam_indx    = BS.beam.data

#Downward looking ADCP
depths       = (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist) - np.cos(beam_angle) * ranges
# Simulated ADCP backscatter data
# Dimensions: time x depth x beam





# Dimensions: time x depth x beam
n_time = len(time)
n_depth = len(depths)
n_beam = len(beam_indx)

#backscatter = np.random.rand(n_time, n_depth, n_beam)  # Random data for demonstration
#depths = np.linspace(0, 100, n_depth)  # Depth levels (0-100 meters)
#times = np.linspace(0, 10, n_time)  # Time points (0-10 seconds)
#beams = np.arange(n_beam)  # Beam index

backscatter = BS  # Random data for demonstration
times       = time  # Time points (0-10 seconds)
beams       = beam_indx  # Beam index

times = mdate.date2num(time)
# Highlighted slice (e.g., at depth index 10)
highlight_depth_idx = 24
highlight_depth = depths[highlight_depth_idx]

#tstart, tend  =  start_endtime(times)
# Create 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Meshgrid for 3D visualization
#T, B = np.meshgrid(times, beams, indexing='ij')
T, B = np.meshgrid(beams, times, indexing='ij')
D = np.full_like(T, highlight_depth)

#print(T)
#print(B)
#exit()
# Plot 3D backscatter slices at each depth
for d in range(n_depth):
    slice_data = backscatter[:, d, :]
    depth_layer = np.full_like(T, depths[d])
    #ax.plot_surface(T, depth_layer, B, rstride=10, cstride=1, facecolors=cm.jet(slice_data / slice_data.max()),
    ax.plot_surface(T, depth_layer, B, rstride=10, cstride=1, facecolors=cm.CMRmap_r(slice_data / slice_data.max()),
                    shade=False, alpha=0.7, edgecolor='none')

# Highlighted slice
highlight_data = backscatter[:, highlight_depth_idx, :]
ax.plot_surface(T, D, B, rstride=10, cstride=1, facecolors=cm.plasma(highlight_data / highlight_data.max()),
#ax.plot_surface(T, D, B, rstride=10, cstride=1, facecolors=cm.CMRmap_r(highlight_data / highlight_data.max()),
                shade=False, alpha=0.9, edgecolor='k')

# Labels and view adjustment
fsize= 13
lpad = 15
lsize = 11
#ax.set_xlabel("Time (s)")
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel("Beam Index", fontsize=fsize, labelpad =lpad)
ax.set_ylabel("Height above lakebed (m)", fontsize=fsize, labelpad =lpad)
#number of vertical lines for grid
locations = verical_grid()
ax.zaxis.set_major_locator(locations)
#########################################3
#Set ylim of x-axis
#ax.set_xlim(tstart, tend)
tstart = datetime.datetime(2023,10,23,0,0)
tend  = datetime.datetime(2023,10,24,0,0)
ax.set_zlim([tstart, tend])
ax.set_zlabel("Time (hour:minute)",fontsize=fsize, labelpad =lpad)
ax.zaxis.set_major_formatter(mdate.DateFormatter("%H:%M"))

##################################################
ax.tick_params(axis ='x', labelsize= lsize)
ax.tick_params(axis ='y', labelsize= lsize)
ax.tick_params(axis ='z', labelsize= lsize)
#ax.set_title("3D ADCP Backscatter with Highlighted Depth Slice")
ax.view_init(elev=30, azim=-45)  # Adjust view angle
#ax.view_init(elev=45, azim=-65)  # Adjust view angle

#plt.colorbar(cm.ScalarMappable(cmap='viridis'), ax=ax, label='Backscatter Intensity')
#plt.colorbar(cm.ScalarMappable(cmap='jet'), ax=ax, label='Backscatter Intensity', pad =0.2, extend ='both', shrink=0.5)
plt.colorbar(cm.ScalarMappable(cmap='CMRmap_r'), ax=ax, label='Backscatter Intensity', pad =0.2, extend ='both', shrink=0.5)
plt.show()
