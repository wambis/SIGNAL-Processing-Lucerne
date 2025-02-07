import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates
import re, pandas as pd 
import dolfyn as dlfn
import datetime
from matplotlib.colors import LinearSegmentedColormap
import json, yaml
from yaml.loader import SafeLoader


#paramerters

Figsize= 14
pad    =10

#Define the grid
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
        depths   = (height_adcp + blank_dist) - np.cos(beam_angle) * ranges
    else:
        depths   = (height_adcp + blank_dist) + np.cos(beam_angle) * ranges
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
            ax.plot_surface(T, B, depth_layer, rstride=10, cstride=1, facecolors=cm.Spectral(slice_data / slice_data.max()),
                            shade=False, alpha=0.7, edgecolor='none')
        
        #Highlighted slice
        highlight_data = backscat[:, depth_idx, :]
        #ax.plot_surface(T, B, D, rstride=10, cstride=1, facecolors=cm.plasma(highlight_data / highlight_data.max()),
        ax.plot_surface(T, B, D, rstride=10, cstride=1, facecolors=cm.PuBu(highlight_data / highlight_data.max()),
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


#Backscatter corrceted function
def BSC_Correct(df, Tlist, ADCP_DOWN= True, MEAN = None):
    #get the parameter
    beam_angle = df.velds.attrs['beam_angle']
    #get the Profile range
    R          = df.range.data
    #Calculate the range along the beam
    #R          = PR * np.cos(beam_angle)
    #########################################
    #make the average on all the 4 beams, NB the key word here is amp
    df_beam_avg = np.mean(df.amp, axis = 0)
    ##############
    #alpha = 0.6
    #############################
    if(MEAN):
        #alpha is not the real value
        alpha = 0.6
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
        #Correct the Bascatter
        PC = P * 0.43 + 20 * np.log10(R) + 2 *alpha * R
    else:
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
        #Get the range, altitude
        #Loop over the list of the time
        Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = R) for ti in Tlist], axis =1 )
        time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
        M2D_NEW     = []
        #Check if the ADCP is downward looking or not
        if(ADCP_DOWN):
            #Frequency 1200 kHz
            alpha = 0.3
            #alpha = 0.48
            for r_i, slide in zip(R, Matrix2D):
                #print(i, slide)
                slide_new = slide *0.43 +  20 * np.log10(r_i) + 2 *alpha * r_i
                #print(slide_new)
                M2D_NEW.append(slide_new)
            #Correct the Bascatter
            PC  = np.asarray(M2D_NEW)
        else:
            #frequency= 600 kHz
            alpha = 0.1
            for r_i, slide in zip(R, Matrix2D):
                #print(i, slide)
                slide_new = slide *0.43 +  20 * np.log10(r_i) + 2 *alpha * r_i
                #print(slide_new)
                M2D_NEW.append(slide_new)
            #Correct the Bascatter
            PC  = np.asarray(M2D_NEW)
    return (time, R, PC)


#def Extract_BS2D(df, beam_indx, Tlist,  param, nrange=None):
#        if(nrange== None):
#            r        = df.range.data
#        else:
#            r        = df.range.data[:nrange]
#        if(param=="backscatter"):
#                #make the average on all the 4 beams, NB the key word here is amp
#                #Get the range, altitude
#                #Loop over the list of the time
#                Matrix2D    = np.concatenate([df.amp.sel(time = ti, range = r) for ti in Tlist], axis =1 )
#                #Matrix2D    = np.concatenate([df.amp.sel(time = ti, range = r) for ti in Tlist])
#                time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
#                return (time, r, Matrix2D[beam_indx])

def Extract_BS2D(df, Tlist,  param, nrange=None):
        if(nrange== None):
            r        = df.range.data
        else:
            r        = df.range.data[:nrange]
        if(param=="backscatter"):
                #make the average on all the 4 beams, NB the key word here is amp
                #Get the range, altitude
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(df.amp, axis = 0)
                #Loop over the list of the time
                #Loop over the list of the time
                #Matrix3D    = np.concatenate([df.amp.sel(time = ti, range = r) for ti in Tlist], axis =1 )
                Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist])
                DList= []
                for i in range(Matrix2D.shape[0]):
                        #print(Matrix3D[beam_indx][i])
                        mean_val = np.mean(Matrix2D[i])
                        List_l = [(item -mean_val)/mean_val for item in Matrix2D[i]]
                        #List_l = [(item -mean_val) for item in Matrix3D[beam_indx][i]]
                        DList.append(List_l)
                #transform DList into an array 
                DList       = np.asarray(DList)
                #print(DList)
                #exit()
                time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
                #return (time, r, Matrix3D[beam_indx])
                return (time, r, DList)


#def Plot_fig2D(ax, fig_object, df, time, data2D,  height_adcp,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
def Plot_fig2D(ax, fig_object, df, time, data2D,  height_adcp,  ADCP_DOWN=False, rframe= None):
    #speed up figure plot
    plt.style.use('fast')
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    #ylabel = "dln(Echo)";  
    ylabel = "Echo Intensity";  
    color = "jet" 
    #color = "RdYlBu_r" 
    # Create custom colormap
    #color = LinearSegmentedColormap.from_list("custom_blue_yellow", colors, N=256)
    #ylabel = 'd'+r'$\ln(Echo)$';  color = "jet" 
    #ylabel = "Backscatter";  color = "jet" 
    vmind         =min(np.min(data2D, axis = 1))
    vmaxd         = max(np.max(data2D, axis = 1))
    ########################
    #vmind = -0.5; vmaxd =+0.5
    #vmind = -0.8; vmaxd =+0.8
    if("velocity" in ylabel):
        vmind = 0.0; vmaxd = 0.4
    #####Test sur L'ADCP orientation
    if(ADCP_DOWN):
        #invert the y-axis
        #ax.invert_yaxis()
        #invert ticks values from positive to negative
        depths    = abs((height_adcp+ blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
        y_lims    = [min(depths), height_adcp]
        #y_lims     = [0, max(depths)]
    else:
        depths    = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
        y_lims    = [height_adcp, max(depths)]
        #y_lims     = [max(depths), 0]
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime'
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    #Set some generic y-limits.
    extent        = [start_num , end_num,  y_lims[0], y_lims[1]]
    #Make a 2D plot
    if(ADCP_DOWN):
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'upper', aspect = 'auto',
                   #resample=True, extent = extent)
                   #interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='bicubic', resample=True, extent = extent)
    else:
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'lower', aspect = 'auto',
                   #resample=True, extent = extent)
                   #interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='bicubic', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #number of vertical lines for grid
    locations   = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks       = ax.get_yticks()
    ##################################
    if('CW' in ylabel):
        #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
        if(ADCP_DOWN):
            #set the color bar on the figure
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.05, cbar_width, cbar_height])
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 90, 180, 270, 350])
            cbar_ax.yaxis.set_label_position("right")
            #shift the label by 3.5 along x-axis and by 1.0 along the y-axis
            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
            #move the label by -0.08 on x-axis and by 1.3 on y-axis
            #ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="top", fontsize = Figsize)
            ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            ax.yaxis.set_label_coords(-0.08 , 1.2)
        #else:
        #    #fig_object.colorbar(im, cax=cbar_ax, ticks =[0, 90, 180, 270, 350])
        #    ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = 12)
    #elif('CW' not in ylabel):
    else:
        if(ADCP_DOWN):
            #set the color bar on the figure
            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.2, cbar_width, cbar_height])
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.1, 0.2, 0.3, 0.4])
            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
            cbar_ax.yaxis.set_label_position("right")
            #cbar_ax.ax.tick_params(labelsize=12)
            cbar_ax.yaxis.set_tick_params(labelsize=12)
            cbar_ax.yaxis.label.set_size(14)
            #position of the color bar
            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
            #ylabel and the postion of its positions
            ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="center", fontsize = Figsize)
            #ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
            #move the label by -0.08 on x-axis and by 1.2 on y-axis
            ax.yaxis.set_label_coords(-0.08 , 1.2)
            #plt.yticks(fontsize = 12)
#        else:
#            #fig_object.colorbar(im, cax=cbar_ax)
#            ax.set_ylabel("H (m)", labelpad = pad, loc="top", fontsize = Figsize)
            #cbar_ax.yaxis.set_label_position("right")
    #Set label
    #Control font size
    #cbar_ax.tick_params(labelsize=10)
    ##############################################################
    #ax.set_yticklabels(ticks)
    #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
    if(ADCP_DOWN):
        #invert the y-axis
        #ax.invert_yaxis()
        ##ax.set_ylim(0, float(height_adcp))
#        ax.invert_yaxis()
        #get yticks
        #yticks = ax.get_yticks()
        ##reverse ticklabels
        #ytick_labels = yticks[::-1]
        #ticks_new = []
        #for ticks in ytick_labels:
        #    if (float(ticks) > height_adcp):
        #        ticks = height_adcp
        #    ticks_new.append(ticks)
        ##ax.set_yticklabels(yticks[::-1])
        ##Set the ticks and labels
        #ax.set_yticklabels(ticks_new)
        ## Optionally, adjust the axis limits to match the ticks
        ##ax.set_ylim(min(ticks_new), max(ticks_new))
        #ax.set_ylim(max(ticks_new), min(ticks_new))
        x=2
    else:
        #ax.set_ylim(float(min(depths)), float(max(depths)))
      #  ax.set_ylim(float(height_adcp), float(max(depths)))
      #  #invert ticks values from positive to negative
      #  ax.annotate(
      #   AphaDict[ix],     # Text for the annotation
      #   #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
      #   xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
      #   xycoords='axes fraction',    # Use axes fraction for positioning
      #   fontsize=16,                 # Font size for the annotation
      #   #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
      #   bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
      #  )
      x=2
    #Set xlim of x-axis
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    #Remove the frame on the plot
    #ax.spines["right"].set_visible(False)
    if(rframe == "top"):
        ax.spines["top"].set_visible(False)
    elif(rframe == "bottom"):
        ax.spines["bottom"].set_visible(False)
        #Remove the ticks marks (small vertical lines) on the x-axis
        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
    elif(rframe == "right"):
        ax.spines["right"].set_visible(False)
    elif(rframe =="left"):
        ax.spines["left"].set_visible(False)
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])









#Open the configuration file
with open("conf2Dbck.yaml") as Fym:
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
#get the date
STARTDATE    = Fig_params['STARTDATE']
ENDDATE      = Fig_params['ENDDATE']
#Grab the dictionary
PARAMSDICT     = Fig_params['PARAMSDICT']
#get the type of slice
VERTICAL_SLIDES  = Fig_params['VERTICAL_SLIDES']
##########################################
ADCP_DOWN      = Fig_params['ADCP_DOWN']
ADCP_UP        = Fig_params['ADCP_UP']
########################################
ADCP_UP_HEIGTH_FROM_LAKEBED     = Fig_params['ADCP_UP_HEIGTH_FROM_LAKEBED']
ADCP_DOWN_HEIGTH_FROM_LAKEBED   = Fig_params['ADCP_DOWN_HEIGTH_FROM_LAKEBED']





try:
    df_down = dlfn.read(ADCP_FILE_NAME_DOWN)
    df_up   = dlfn.read(ADCP_FILE_NAME_UP)
except:
    print("print provide the ADCP files:  %s and %s  "%(ADCP_FILE_NAME_UP, ADCP_FILE_NAME_DOWN))


try:
    date_all = pd.date_range(start= STARTDATE, end = ENDDATE)
except:
     print("Checked the date entring %s  %s"%(STARTDATE, ENDDATE))
     exit()
#Create a time list for the started and the
date_list = [str(it).split()[0] for it in date_all]



#get the backscatter
if(ADCP_UP==True and ADCP_DOWN==True):
    nfigs = 2
elif(ADCP_UP==True and ADCP_DOWN==False or ADCP_UP==False and ADCP_DOWN==True):
    nfigs = 1


#Grab the space between the figure
fig_space = Fig_params['fig_space']

#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

#Create the figure
fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)




if(ADCP_UP==True and ADCP_DOWN==True):
        param = 'backscatter'
        #for bindx in range(4):
        time, rup, Matrix2D_up      = BSC_Correct(df_up, date_list, ADCP_DOWN= False)
        time, rdown, Matrix2D_down  = BSC_Correct(df_down, date_list, ADCP_DOWN= True)
        #plot the 2D Backscatter
        #plot the upward looking ADCP
        Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #plot the downward looking ADCP
        Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
elif(ADCP_UP==True and ADCP_DOWN==False):
        param = 'backscatter'
        time, r, Matrix2D  = BSC_Correct(df_up, date_list, ADCP_DOWN= False)
        Plot_fig2D(axs, fig, df_up, time, Matrix2D, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)

elif(ADCP_UP==False and ADCP_DOWN==True):
        param   = 'backscatter'
        time, r, Matrix2D   = BSC_Correct(df_down, date_list, ADCP_DOWN= True)
        #plot the downward looking ADCP
        Plot_fig2D(axs, fig, df_down, time, Matrix2D, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        print(Matrix2D.shape)
#################### start ###########################
plt.subplots_adjust(hspace = fig_space)
#######
if(nfigs>1):
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
else:
    axs.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

figname = "Figure_Backscatter2D_Corrected%s.png"%(STARTDATE)
#Save the figure
plt.savefig(figname, bbox_inches = 'tight')
