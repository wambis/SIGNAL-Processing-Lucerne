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
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
from yaml.loader import SafeLoader
######Oceanic colors
import cmocean

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




def Plot_Temp2D(ax, fig_object, time, data2D, depths, **PARAMSDICT):
    #speed up figure plot
    plt.style.use('fast')
    ##Get the parameters
    ylabel, color, ix = PARAMSDICT[param]
    vmind         = np.min(data2D, axis = 1)
    vmaxd         = np.max(data2D, axis = 1)
    #get the start and the endtime
    tstart, tend  =  start_endtime(time)
    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime'
    startdatetime = pd.to_datetime(tstart).to_pydatetime()
    enddatetime   = pd.to_datetime(tend).to_pydatetime()
    ###########################################
    start_num     = mdates.date2num(startdatetime)
    end_num       = mdates.date2num(enddatetime)
    #Set some generic y-limits.
    y_lim_min     = np.min(depths)
    y_lim_max     = np.max(depths)
    #Set the extent for 2D plot
    extent        = [start_num , end_num,  y_lim_min, y_lim_max]
    #print(extent)
    #Make a 2D plot
    im            = ax.imshow(data2D, cmap = color, vmin = min(vmind) , vmax = max(vmaxd), origin = 'lower', aspect = 'auto',
                   interpolation ='spline16', resample=True, extent = extent)
                   #interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
    #get the minimum and the maximum of yaxis
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    #get the coordinate of the point, where to write the text
    xp, yp     = xy_point(xmin, xmax,ymin, ymax)
    #Write the text on the figure
#    ax.annotate(AphaDict[ix] , xy=(xp, yp), fontsize = f_size)
    ax.annotate(
         AphaDict[ix],     # Text for the annotation
         #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
         xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
         xycoords='axes fraction',    # Use axes fraction for positioning
         fontsize=16,                 # Font size for the annotation
         #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
         bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
    )
    #Add the colorbar using the figure's method,telling which mappable we're talking about and
    #which axes object it should be near Calculate the colorbar position and size
    bbox        = ax.get_position()
    cbar_width  = 0.02
    cbar_height = bbox.height
    cbar_x      = bbox.x1 + 0.02
    cbar_y      = bbox.y0
    #Create Colorbar's axis
    cbar_ax  = fig_object.add_axes([cbar_x , cbar_y, cbar_width, cbar_height])
    ########################
    #Add the colorbar on the figure
    fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
    #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[6, 6.5,7 ,7.5,8.5,9])
    #################################
    #Set label
    ax.set_ylabel('H (m)', labelpad = pad, fontsize = Figsize)
    #into a nice datetime string.
    #ax.xaxis_date()
    #number of vertical lines for grid
    locations = verical_grid()
    ax.xaxis.set_major_locator(locations)
    #get ticks values
    ticks         = ax.get_yticks()
    #ax.set_yticklabels(ticks)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(base=20))
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    ax.set_ylim(float(min(depths)), float(max(depths)))
    #Set xlim of x-axis
    ax.set_xlim(tstart, tend)
    #Make the grid
    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #disable ticks on  x-axis
    ax.set_xticklabels([])






#def Scatt_Correct(df, Profile_range, angle = 20.0):
def BSC_Correct(df, Tlist, alpha =0.6, nrange=25, MEAN = True):
    #get the parameter
    beam_angle = df.velds.attrs['beam_angle']
    #get the Profile range
    PR         = df.range.data[: nrange]
    #Calculate the range along the beam
    R          = PR * np.cos(beam_angle)
    #########################################
    #amp_new    = df.amp * 0.43 + 20 * np.log10(R) + 2 *alpha * R
    #alpha is water absorption
    #alpha      = 0.6
    #alpha      = 2.4  #Good value for both ADCP (UP and down)
    ####################################
    #make the average on all the 4 beams, NB the key word here is amp
    df_beam_avg = np.mean(df.amp, axis = 0)
    #############################
    if(MEAN):
        #Loop over the list of the time
        P  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in Tlist])
        T  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in Tlist])
        #Correct the Bascatter
        PC = P * 0.43 + 20 * np.log10(R) + 2 *alpha * R
    else:
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
        #Get the range, altitude
        #r           =  df_beam_avg.range
        r           =  df_beam_avg.range.data[:nrange]
        #Loop over the list of the time
        Matrix2D    = np.concatenate([df_beam_avg.sel(time = ti, range = r) for ti in Tlist], axis =1 )
        #time        = np.concatenate([df.amp.sel(time = ti) for ti in Tlist])
        time        = np.concatenate([df.amp.sel(time = ti)['time'] for ti in Tlist])
        M2D_NEW     = []
        #for ir, slide in zip(np.arange(len(Matrix2D)), Matrix2D):
        for r_i, slide in zip(R, Matrix2D):
            #print(i, slide)
            slide_new = slide *0.43 +  20 * np.log10(r_i) + 2 *alpha * r_i
            #print(slide_new)
            M2D_NEW.append(slide_new)
            #print(type(slide_new))
        #Correct the Bascatter
        #PC = Matrix2D * 0.43 + 20 * np.log10(R) + 2 *alpha * R
        PC  = np.asarray(M2D_NEW)
    return time, PC


###############################################
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
                #mean_val = np.mean(np.mean(df_beam_avg, axis = 1), axis = 0)
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

#
##def Plot_fig2D(ax, fig_object, df, time, data2D,  height_adcp,  ADCP_DOWN=False, rframe= None, **PARAMSDICT):
#def Plot_fig2D(ax, fig_object, df, time, data2D,  height_adcp,  ADCP_DOWN=False, rframe= None):
#    #speed up figure plot
#    plt.style.use('fast')
#    data2D = data2D * 100
#    #get the parameters
#    beam_angle   = df.beam_angle
#    blank_dist   = df.blank_dist
#    ranges       = df.range.data
#    ##Get the parameters
#    ylabel = "dln(Echo)";  
#    #color = "jet" 
#    #color = plt.get_cmap('autumn_r') 
#    #color = cmocean.cm.delta
#    #color = cmocean.cm.thermal
#    #color = cmocean.cm.haline
#    #color = cmocean.cm.deep_r
#    #color = cmocean.cm.speed_r
#    #color = cmocean.cm.rain_r
#    #color = cmocean.cm.curl
#    color = cmocean.cm.tarn_r
#    #color = "RdYlBu_r" 
#    # Create custom colormap
#    #color = LinearSegmentedColormap.from_list("custom_blue_yellow", colors, N=256)
#    #ylabel = 'd'+r'$\ln(Echo)$';  color = "jet" 
#    #ylabel = "Backscatter";  color = "jet" 
#    vmind         =min(np.min(data2D, axis = 1))
#    vmaxd         = max(np.max(data2D, axis = 1))
#    ########################
#    vmind = -50; vmaxd =+50
##    vmind = -50; vmaxd =+70
#    #vmind = -70; vmaxd =+60
#    #vmind = -0.5; vmaxd =+0.5
#    #vmind = -0.8; vmaxd =+0.8
#    if("velocity" in ylabel):
#        vmind = 0.0; vmaxd = 0.4
#    #####Test sur L'ADCP orientation
#    if(ADCP_DOWN):
#        #invert the y-axis
#        #ax.invert_yaxis()
#        #invert ticks values from positive to negative
#        depths    = abs((height_adcp+ blank_dist) - np.cos(np.radians(beam_angle)) * ranges)
#        y_lims    = [min(depths), height_adcp]
#        #y_lims     = [0, max(depths)]
#    else:
#        depths    = (height_adcp + blank_dist) + np.cos(np.radians(beam_angle)) * ranges
#        y_lims    = [height_adcp, max(depths)]
#        #y_lims     = [max(depths), 0]
#        #Remove the ticks marks (small vertical lines) on the x-axis
#        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
#    #get the start and the endtime
#    tstart, tend  =  start_endtime(time)
#    ###### Convert 'numpy.datetime64' ==>  'datetime.datetime'
#    startdatetime = pd.to_datetime(tstart).to_pydatetime()
#    enddatetime   = pd.to_datetime(tend).to_pydatetime()
#    ###########################################
#    start_num     = mdates.date2num(startdatetime)
#    end_num       = mdates.date2num(enddatetime)
#    #Set some generic y-limits.
#    extent        = [start_num , end_num,  y_lims[0], y_lims[1]]
#    #Make a 2D plot
#    if(ADCP_DOWN):
#        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'lower', aspect = 'auto',
#                   #interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
#                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
#                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)
#    else:
#        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'upper', aspect = 'auto',
#                   #interpolation ='spline16',interpolation_stage='rgba', resample=True, extent = extent)
#                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
#                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)
#        cmocean.plots.test(cmocean.cm.gray, ax=ax)
#    #get the minimum and the maximum of yaxis
#    ymin, ymax = ax.get_ylim()
#    xmin, xmax = ax.get_xlim()
#    #get the coordinate of the point, where to write the text
#    #xp, yp     = xy_point(xmin, xmax,ymin, ymax)
#    #which axes object it should be near Calculate the colorbar position and size
#    bbox        = ax.get_position()
#    cbar_width  = 0.02
#    cbar_height = bbox.height
#    cbar_x      = bbox.x1 + 0.02
#    cbar_y      = bbox.y0
#    #number of vertical lines for grid
#    locations   = verical_grid()
#    ax.xaxis.set_major_locator(locations)
#    #get ticks values
#    ticks       = ax.get_yticks()
#    ##################################
#    if(ADCP_DOWN):
#            #set the color bar on the figure
#            cbar_ax  = fig_object.add_axes([cbar_x , cbar_y +0.2, cbar_width, cbar_height])
#            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
#            #fig_object.colorbar(im, cax=cbar_ax, label = ylabel,ticks =[0, 0.1, 0.2, 0.3, 0.4])
#            fig_object.colorbar(im, cax=cbar_ax, label = ylabel)
#            cbar_ax.yaxis.set_label_position("right")
#            #cbar_ax.ax.tick_params(labelsize=12)
#            cbar_ax.yaxis.set_tick_params(labelsize=12)
#            cbar_ax.yaxis.label.set_size(14)
#            #position of the color bar
#            #cbar_ax.yaxis.set_label_coords(3.5, 1.08)
#            #ylabel and the postion of its positions
#            ax.set_ylabel("Height above lakebed (m)", labelpad = pad, loc="center", fontsize = Figsize)
#            #ax.set_ylabel("H (m)", labelpad = pad, loc="center", fontsize = Figsize)
#            #move the label by -0.08 on x-axis and by 1.2 on y-axis
#            ax.yaxis.set_label_coords(-0.08 , 1.2)
#            #plt.yticks(fontsize = 12)
##        else:
##            #fig_object.colorbar(im, cax=cbar_ax)
##            ax.set_ylabel("H (m)", labelpad = pad, loc="top", fontsize = Figsize)
#            #cbar_ax.yaxis.set_label_position("right")
#    #Set label
#    #Control font size
#    #cbar_ax.tick_params(labelsize=10)
#    ##############################################################
#    #ax.set_yticklabels(ticks)
#    #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
#    if(ADCP_DOWN):
#        #invert the y-axis
#        #ax.invert_yaxis()
#        ##ax.set_ylim(0, float(height_adcp))
##        ax.invert_yaxis()
#        #get yticks
#        #yticks = ax.get_yticks()
#        ##reverse ticklabels
#        #ytick_labels = yticks[::-1]
#        #ticks_new = []
#        #for ticks in ytick_labels:
#        #    if (float(ticks) > height_adcp):
#        #        ticks = height_adcp
#        #    ticks_new.append(ticks)
#        ##ax.set_yticklabels(yticks[::-1])
#        ##Set the ticks and labels
#        #ax.set_yticklabels(ticks_new)
#        ## Optionally, adjust the axis limits to match the ticks
#        ##ax.set_ylim(min(ticks_new), max(ticks_new))
#        #ax.set_ylim(max(ticks_new), min(ticks_new))
#        x=2
#    else:
#        #ax.set_ylim(float(min(depths)), float(max(depths)))
#      #  ax.set_ylim(float(height_adcp), float(max(depths)))
#      #  #invert ticks values from positive to negative
#      #  ax.annotate(
#      #   AphaDict[ix],     # Text for the annotation
#      #   #xy=(0.0089, 0.903),                   # Coordinates (relative to axes) for the annotation
#      #   xy=(0.0089, 0.72),                   # Coordinates (relative to axes) for the annotation
#      #   xycoords='axes fraction',    # Use axes fraction for positioning
#      #   fontsize=16,                 # Font size for the annotation
#      #   #bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5)  # Box style
#      #   bbox=dict(boxstyle='square', facecolor='white', edgecolor='k', alpha=0.5)  # Box style
#      #  )
#      x=2
#    #Set xlim of x-axis
#    #Make the grid of the axis
#    ax.grid(visible = True, axis = "x", alpha = 0.2)
#    # Make ticks on x-axis and y-axis bigger
#    #plt.tick_params(axis='both', which='major', labelsize=14)
#    ax.tick_params(axis='both', which='major', labelsize=14)
#    #Remove the frame on the plot
#    #ax.spines["right"].set_visible(False)
#    if(rframe == "top"):
#        ax.spines["top"].set_visible(False)
#    elif(rframe == "bottom"):
#        ax.spines["bottom"].set_visible(False)
#        #Remove the ticks marks (small vertical lines) on the x-axis
#        ax.tick_params(axis="x", length = 0, color= "white", width = 0)
#    elif(rframe == "right"):
#        ax.spines["right"].set_visible(False)
#    elif(rframe =="left"):
#        ax.spines["left"].set_visible(False)
#    #Set xlim of x-axis
#    ax.set_xlim(tstart, tend)
#    #Make the grid
#    #ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#    #disable ticks on  x-axis
#    ax.set_xticklabels([])


def Extract_RBR_SOLO_TEMPERATURE(DataDict,  Tlist, AVERAGE=False):
    #The data should be in the dictionary, where the keys are the depths
    DATA2D     = []
    DEPTHS     = []
    #get the time you need to plot
    TIMES      = set([time for depth in DataDict.keys() for time in DataDict[depth] if(time.split()[0] in Tlist)])
    #####################################
    TIMES      = sorted(np.asarray(list(TIMES)))
    #loop over the depths and the time to form a matrix
    #depthss = [d for d in DataDict.keys()]
    #depths_up   = [70, 60, 50, 37.5]
    #depths_down = [20, 17.5, 15, 5 ]

    depths_up   = [40,  55]
    depths_down = [5, 12.5]
    ##################################
#    depths_up   = [70]
#    depths_down = [20]
    #['70', '67.5', '65', '62.5', '60', '57.5', '55', '52.5', '50', '47.5', '45', '42.5', '40', '37.5', '35', '32.5', '30', '27.5', '25', '22.8', '20', '17.5',  
    #'15', '12.5', '10', '7.5', '5']
    #print(depthss)
    DICT_UP   = dict()
    DICT_DOWN = dict()
    #exit()
    #for  depth in DataDict.keys():
    for  depth in depths_up:
        #get the temperature at a specific depth
        temp   = [DataDict[str(depth)][time] for time in TIMES]
        #####################################
        DICT_UP[depth] = (pd.to_datetime(TIMES), temp)
    ########### loop over the depth down
    for  depth in depths_down:
        #get the temperature at a specific depth
        temp   = [DataDict[str(depth)][time] for time in TIMES]
        #####################################
        DICT_DOWN[depth] =(pd.to_datetime(TIMES), temp)
    #Trans formed the List into# 2D array
    return(DICT_UP, DICT_DOWN)




def Extract_RBR_SOLO_TEMP(DataDict,  Tlist, AVERAGE=False):
    #The data should be in the dictionary, where the keys are the depths
    DATA2D     = []
    DEPTHS     = []
    #get the time you need to plot
    TIMES      = set([time for depth in DataDict.keys() for time in DataDict[depth] if(time.split()[0] in Tlist)])
    #####################################
    TIMES      = sorted(np.asarray(list(TIMES)))
    #loop over the depths and the time to form a matrix
    for  depth in DataDict.keys():
        #get the temperature at a specific depth
        temp   = [DataDict[depth][time] for time in TIMES]
        #print(temp)
        DATA2D.append(temp)
        DEPTHS.append(float(depth))
    #Trans formed the List into# 2D array
    DATA2D     = np.asarray(DATA2D, dtype=float)
    #Convert the TIMES into pan#das time
    TIMES      = pd.to_datetime(TIMES)
    #Check if you need everage #at the depth
    if(AVERAGE):
        data   = np.mean(DATA2D, axis = 0)
        return (TIMES, data)
    else:
        return(TIMES,DEPTHS, DATA2D)










#def Plot_fig2D(ax, fig_object, df, time, data2D, temp_data2D, height_adcp,  ADCP_DOWN=False, rframe= None):
def Plot_fig2D(ax, fig_object, df, time, data2D, Dictdata2D, height_adcp,  ADCP_DOWN=False, rframe= None):
    #speed up figure plot
    plt.style.use('fast')
    data2D = data2D * 100
    #get the parameters
    beam_angle   = df.beam_angle
    blank_dist   = df.blank_dist
    ranges       = df.range.data
    ##Get the parameters
    ylabel = "dln(Echo)";  
    #color = "jet" 
    #color = cmocean.cm.thermal
    color = cmocean.cm.haline
    #color = cmocean.cm.deep_r
    #color = cmocean.cm.speed_r
    #color = cmocean.cm.rain_r
    #color = cmocean.cm.curl
    #color = "RdYlBu_r" 
    ############################
    # Create custom colormap
    #color = LinearSegmentedColormap.from_list("custom_blue_yellow", colors, N=256)
    #ylabel = 'd'+r'$\ln(Echo)$';  color = "jet" 
    #ylabel = "Backscatter";  color = "jet" 
    vmind         =min(np.min(data2D, axis = 1))
    vmaxd         = max(np.max(data2D, axis = 1))
    ########################
    vmind = -50; vmaxd =+50
#    vmind = -50; vmaxd =+70
    #vmind = -70; vmaxd =+60
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
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'lower', aspect = 'auto',
                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)
        #####Temperature Contour##############
        #contour = ax.contour(time, data2D, levels=10, colors='black',extent= extent )
        #contour       = ax.contourf(times, periods, np.abs(coefs), levels=10, cmap='viridis', extent= extent)

        #ax2 = ax.twinx()
        #for temp in temp_data2D:
        #    #temp = np.asarray(Dictdata2D[key][1])
        #    ax2.plot(time, temp, color='k', lw=2.5)
        #    #ax2.set_ylim(depth_to_plot, depth_to_plot)  # Set the y-axis of ax2 to just the specific depth
        #    #ax2.set_ylim(70, 70)  # Set the y-axis of ax2 to just the specific depth
                   #
    else:
        im        = ax.imshow(data2D, cmap = color, vmin = vmind , vmax = vmaxd, origin = 'upper', aspect = 'auto',
                   #interpolation ='bicubic',interpolation_stage='rgba', resample=True, extent = extent)
                   interpolation ='hanning',interpolation_stage='rgba', resample=True, extent = extent)
        #####Temperature Contour##############
        #for temp in temp_data2D:
            #temp = np.asarray(Dictdata2D[key][1])
        #    ax.plot(time, temp, color='k', lw=2.5)
        #contour = ax.contour(XX, YY, data2D, levels=10, colors='black')
        #ax.clabel(contour, inline=True, fontsize=18)
    #get the minimum and the maximum of yaxis
    ymin, ymax  = ax.get_ylim()
    xmin, xmax  = ax.get_xlim()
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
    # Add a label or annotation at a specific time, e.g., at 01:00
    target_time = datetime(2023, 10, 23, 2, 24)  # 02:24
    #cbar_ax.tick_params(labelsize=10)
    ax2 = ax.twinx()  # Create a second axis sharing the same x-axis
    for di in Dictdata2D.keys():
            timei     = Dictdata2D[di][0]
            tempMean  = np.mean(np.asarray(Dictdata2D[di][1]))
            temp      = np.asarray(Dictdata2D[di][1]) + float(di)
            ax2.plot(timei, temp, color='white', lw=1.0)
            #ax2.axhline(float(di), color='r', linestyle='--', label='%f reference'%(di))
            ax.text(target_time,  di+ tempMean *1.1,  "%.2f °C"%(tempMean), color="white", fontsize=12)
           # print(di, tempMean)
    #remove label on y-ax2 axis
    #ax2.set_yticks([])  # Removes ticks
    #ax2.set_yticklabels([])  # Removes the labels
    ax2.axis('off')
    ##############################################################
    # Set tick spacing to 10
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax2.yaxis.set_major_locator(MultipleLocator(5))
    #ax.set_yticklabels(ticks)
    #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_reverse))
    if(ADCP_DOWN):
        #ax.set_ylim(-0.1, float(height_adcp))
        #ax2.set_ylim(-0.1, float(height_adcp))

        ax.set_ylim(float(min(depths)), float(max(depths)))
        ax2.set_ylim(float(min(depths)), float(max(depths)))
        #get yticks
    else:
        ax.set_ylim(float(min(depths)), float(max(depths)))
        ax2.set_ylim(float(min(depths)), float(max(depths)))
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
    #Set xlim of x-axis
    #Make the grid of the axis
    ax.grid(visible = True, axis = "x", alpha = 0.2)
    # Make ticks on x-axis and y-axis bigger
    #plt.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='major', labelsize=14)
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

#grab the temperature file
FILE_TEMPE     = Fig_params['FILE_TEMPE']

#Grab the space between the figure
fig_space = Fig_params['fig_space']

#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])

#Create the figure
fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)



#Grab 2D Temperature
#get the temperature data from RBR-SOLO-3
# Opening JSON file
ftmp   = open(FILE_TEMPE)
# returns JSON object as
# a dictionary
data_temp = json.load(ftmp)
#get the horizontal velocity data for the corresponding time
#time, depths, temp_data2D       = Extract_RBR_SOLO_TEMP(data_temp,  date_list, AVERAGE=False)
#time, depths, temp_data2D       = Extract_RBR_SOLO_TEMP(data_temp,  date_list, AVERAGE=False)
DICT_UP, DICT_DOWN   = Extract_RBR_SOLO_TEMPERATURE(data_temp,  date_list, AVERAGE=False)

depths_up   = [70, 60, 50, 37.5]
#depths_up   = 70
depths_down = [20, 17.5, 15, 5 ]
###########
#index_up = []
#index_down = []
##Find the index of the depth of interest
#for depth in depths_up:
#    depth_index = np.abs(np.asarray(depths) - depth).argmin()
#    index_up.append(depth_index)
#
#for depth in depths_down:
#    depth_index = np.abs(np.asarray(depths) - depth).argmin()
#    index_down.append(depth_index)
#
##Extract the temperature at the specified depth
#temps_up   = temp_data2D[index_up, :]
#temps_down = temp_data2D[index_down, :]
##print(depth_index)

#get the parameter, the index of the corresponding axis, and the color
#_ , color, ix            = PARAMSDICT[param]
#Plot the figure by calling the plotting function of 2D Matrix
#Plot_Temp2D(axs[ix],fig, time, data2D, depths,  **PARAMSDICT)


if(ADCP_UP==True and ADCP_DOWN==True):
        param = 'backscatter'
        #for bindx in range(4):
        time, rup, Matrix2D_up     = Extract_BS2D(df_up,  date_list,  param, nrange=None)
        time, rdown, Matrix2D_down = Extract_BS2D(df_down,  date_list,  param, nrange=None)
        #plot the 2D Backscatter
        #plot the upward looking ADCP
        #Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, temps_up, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        Plot_fig2D(axs[0], fig, df_up, time, Matrix2D_up, DICT_UP, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)
        #plot the downward looking ADCP
        #Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        #Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, temps_down, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
        Plot_fig2D(axs[1], fig, df_down, time, Matrix2D_down, DICT_DOWN, ADCP_DOWN_HEIGTH_FROM_LAKEBED, ADCP_DOWN=True, rframe= None)
elif(ADCP_UP==True and ADCP_DOWN==False):
        param = 'backscatter'
        time, r, Matrix2D  = Extract_BS2D(df_up,  date_list,  param)
        Plot_fig2D(axs, fig, df_up, time, Matrix2D, ADCP_UP_HEIGTH_FROM_LAKEBED, ADCP_DOWN=False, rframe= None)

elif(ADCP_UP==False and ADCP_DOWN==True):
        param   = 'backscatter'
        time, r, Matrix2D   = Extract_BS2D(df_down, date_list,  param)
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

figname = "Figure_Backscatter2D_%s.png"%(STARTDATE)
#Save the figure
plt.savefig(figname, bbox_inches = 'tight')
