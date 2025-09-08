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
import json, yaml
from yaml.loader import SafeLoader


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



#def Plot3DBackscatter(ax, df,TARGETDATE, height_adcp, depth_idx, ADCP_DOWN=False, **PARAMSDICT):
def Plot3DBackscatter(ax, df,TARGETDATE, height_adcp, depth_idx, ADCP_DOWN=False, VERTICAL_SLIDES =None):
    #get the backscatter for the target date
    backscat     = df.amp.sel(time = TARGETDATE)
    #get the paramter
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = backscat.range.data
    time         = backscat.time.data
    beam_indx    = backscat.beam.data
    #Check the if you're plotting the upward, downlard
    if(ADCP_DOWN):
        depths   = abs( (height_adcp + blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
    else:
        depths   = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
    #Dimensions: time x depth x beam
    n_time  = len(time)
    n_depth = len(depths)
    n_beam  = len(beam_indx)
    ################################# 
    beams   = beam_indx  # Beam index
    #get the time in the matplotlib type 
    times = mdate.date2num(time)
    # Highlighted slice (e.g., at depth index 10)
    highlight_depth = depths[depth_idx]
    #Plot the verticals slice
    if(VERTICAL_SLIDES):
        T, B = np.meshgrid(beams, times, indexing='ij')
        D    = np.full_like(T, highlight_depth)
        #loop over the depths
        for d in range(n_depth):
            slice_data  = backscat[:, d, :]
            depth_layer = np.full_like(T, depths[d])
            #ax.plot_surface(T, B, depth_layer, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
            ax.plot_surface(B,T, depth_layer, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
                            shade=False, alpha=0.7, edgecolor='none')
        
        #Highlighted slice
        highlight_data = backscat[:, depth_idx, :]
        #ax.plot_surface(T, B, D, rstride=10, cstride=1, facecolors=cm.PuBu(highlight_data / highlight_data.max()),
        ax.plot_surface(B, T, D, rstride=10, cstride=1, facecolors=cm.PuBu(highlight_data / highlight_data.max()),
                        shade=False, alpha=0.9, edgecolor='k')
        
    #Plot the horizontal slice
    else:
            T, B = np.meshgrid(beams, times, indexing='ij')
            D    = np.full_like(T, highlight_depth)
            #Plot 3D backscatter slices at each depth
            for d in range(n_depth):
                slice_data = backscat[:, d, :]
                depth_layer = np.full_like(T, depths[d])
                ax.plot_surface(T, depth_layer, B, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
                                shade=False, alpha=0.7, edgecolor='none')
            #Highlighted slice
            highlight_data = backscat[:, depth_idx, :]
            ax.plot_surface(T, D, B, rstride=10, cstride=1, facecolors=cm.plasma(highlight_data / highlight_data.max()),
                            shade=False, alpha=0.9, edgecolor='k')




#Open the configuration file
with open("confbck.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)

#Grab the path of the data
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
#Discharge file name
#Grab the Meteo data file
#############################################
#Check the user want to make a plot bigger?
plot_bigger    = Fig_params['plot_bigger']
#Grab the starttime
TARGETDATE     = Fig_params['TARGETDATE']
#Grab the dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#get the type of slice
VERTICAL_SLIDES  = Fig_params['VERTICAL_SLIDES']
##########################################
ADCP_DOWN      = Fig_params['ADCP_DOWN']
#ADCP_UP        = Fig_params['ADCP_UP']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED     = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']

#get the depth index the user want to plot
depth_idx    = Fig_params['depth_idx'] 


if(ADCP_DOWN):
    df           = dlfn.read(ADCP_FILE_NAME_DOWN) 
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    depths       = abs( (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
elif(ADCP_UP):
    df           = dlfn.read(ADCP_FILE_NAME_UP) 
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    depths       = (ADCP_DOWN_HEIGTH_FROM_LAKEBED + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
else:
    print("print provide the ADCP file %s "%(ADCP_FILE_NAME_DOWN))
#get the backscatter


#Create 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
#Labels and view adjustment
fsize= 13
lpad = 12
lsize = 11
#set this for time limit
tstart = datetime.datetime(2023,10,23,0,0)
tend   = datetime.datetime(2023,10,24,0,0)
#Check if the user want to plot the vertical slice
if(VERTICAL_SLIDES):
    figname = "Backscatter_3D_Vertical_Slices.png"
    #figname = "Backscatter_3D_Vertical_Slices.pdf"
    if(ADCP_DOWN):
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_DOWN_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    else:
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_UP_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_ylabel("Beams", fontsize=fsize, rotation=100,labelpad =lpad)
    ax.set_ylabel("Beams", fontsize=fsize)
    ax.set_zlabel("Height above lakebed (m)", fontsize=fsize, labelpad =lpad)
    ax.set_zlim(top=max(depths))
    #ax.set_zlim(min(depths), max(depths))
    #number of vertical lines for grid
    locations    = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #control the visibility of the grid
    ax.set_box_aspect([1,1,1])
    ax.grid(False)
    ax.set_facecolor((1,1,1,0))
    ax.xaxis.pane.fill = False  #Remove pane background
    ax.yaxis.pane.fill = False  #Remove pane background
    ax.zaxis.pane.fill = False  #Remove pane background
    ########################################3
    #Set ylim of x-axis
    ax.set_xlim([tstart, tend])
    #ax.set_ylabel("Time (hour:minute)",fontsize=fsize, labelpad =lpad)
    ax.set_xlabel("Time (hour:minute)",fontsize=10, labelpad =12)
    ax.xaxis.set_major_formatter(mdate.DateFormatter("%H:%M"))
    ##################################################
    #ax.tick_params(axis ='y', labelsize= lsize, pad=1.5)
    ax.tick_params(axis ='y', labelsize= lsize, pad=1.5)
    ax.tick_params(axis ='x', labelsize= 9, labelrotation=310,zorder=12, pad=0.5)
    ax.tick_params(axis ='z', labelsize= lsize,zorder=12)
    #ax.set_title("3D ADCP Backscatter with Highlighted Depth Slice")
    #ax.view_init(elev=20, azim=-30)  # Adjust view angle
    ax.view_init(elev=9.5, azim=-82)  # Adjust view angle
else:
    figname = "Backscatter_3D_Horizontal_Slices.png"
    pad_ticks_h  = 1.2
    if(ADCP_DOWN):
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_DOWN_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    else:
        Plot3DBackscatter(ax, df,TARGETDATE, ADCP_UP_HEIGTH_FROM_LAKEBED, depth_idx, ADCP_DOWN, VERTICAL_SLIDES)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_xlabel("Beam index", fontsize=fsize, labelpad =lpad)
    ax.set_xlabel("Beams", fontsize=fsize, labelpad =lpad)
    ax.set_ylabel("Height above lakebed (m)", fontsize=fsize, labelpad =lpad)
    #number of vertical lines for grid
    locations = verical_grid()
    ax.zaxis.set_major_locator(locations)
    #########################################3
    #Set zlim of x-axis
    ax.set_zlim([tstart, tend])
    ax.set_zlabel("Time (hour:minute)",fontsize=fsize, labelpad =lpad)
    ax.zaxis.set_major_formatter(mdate.DateFormatter("%H:%M"))
    ##################################################
    ax.tick_params(axis ='x', labelsize= lsize, pad= pad_ticks_h)
    ax.tick_params(axis ='y', labelsize= lsize, pad= pad_ticks_h)
    ax.tick_params(axis ='z', labelsize= lsize)
    #ax.set_title("3D ADCP Backscatter with Highlighted Depth Slice", loc="center", pad = -15.0)
    #ax.view_init(elev=30, azim=-45)  # Adjust view angle
    ####################
    #ax.view_init(elev=20, azim=-30)  # Adjust view angle
    ax.view_init(elev=20, azim=10)  # Adjust view angle
#Add the colorbar and save the figure
plt.colorbar(cm.ScalarMappable(cmap='Spectral'), ax=ax, label='Backscatter Intensity', pad =0.04, extend ='both', shrink=0.5)
plt.savefig(figname, bbox_inches = 'tight')
