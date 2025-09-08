import pandas as pd
import numpy as np
from windrose import WindroseAxes
import matplotlib.pyplot as plt
import os, gc, json, yaml
from yaml.loader import SafeLoader
#################################
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import xarray as xr
import pandas as pd
from mhkit import dolfyn as dlfn
#################
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
import matplotlib.cm as cm
import matplotlib.colors as mcolors

#Figsize= 18
Figsize = 16
pad     = 20
#Anotation fontsize
FsizeA  = 20
###############
#############################
def check_if_string(var):
    return isinstance(var, str)
#################################
def inter1DD(x, y, x_new):
    #Compute the interpolat
    ff   = interpolate.interp1d(x, y)
    #define the new points
    xnew = np.linspace(min(x), max(x), num=len(x_new))
    ynew = ff(xnew)   # use interpolation function returned by `interp1d`
    #return the new values
    return(xnew, ynew)
#################

def Extract_ADCPData(df,  Tlist,  param, average = None):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    #Create an empty list to append the extracted data Check if the parameter==velocity
    #Create an empty list to append the extracted data Check if the parameter==velocity
    if(check_if_string(Tlist)):
        Tlist =[ Tlist ]
    #Check if the desired velocity is the Horizontal velocity
    if(param=="velocity_up" or param=="velocity_down"):
        #######################
        #Loop over the list of the time
        try:
            #Loop over the list of the time
            P_list = [np.nanmean(df.velds.U_mag.sel(time=ti), axis=0) for ti in Tlist]
            T_list = [df.velds.U_mag.sel(time=ti)['time'].values for ti in Tlist]
            ###############
            # Flatten lists and create DataFrame
            df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
            # Averaging duplicate timestamps
            df_extracted = df_extracted.groupby('time', as_index=False).mean()
            ####################################
            P            = df_extracted['value']
            T            = df_extracted['time'] 
            #################################
            del P_list ,T_list ,df_extracted 
            #clean the memory
            gc.collect()  # Force garbage collection
        except:
            print("The Data seems not to exist for the certain data in the list Time :  ")
            print(' , '.join(Tlist))
            exit()

    #Check if the desired velocity is the vertical velocity
    elif(param=="vertical_vel"):
        #Loop over the list of the time
        #Vertical component Z--vertical componet
#        P_list = [np.nanmean(df.velds.w.sel(time=ti), axis=0) for ti in Tlist]
#        T_list = [df.velds.w.sel(time=ti)['time'].values for ti in Tlist]
#        #North-Sud component N-S Component
        P_list = [np.nanmean(df.velds.v.sel(time=ti), axis=0) for ti in Tlist]
        T_list = [df.velds.v.sel(time=ti)['time'].values for ti in Tlist]
        #West-Est component  W-E Component
        #P_list = [np.nanmean(df.velds.u.sel(time=ti), axis=0) for ti in Tlist]
        #T_list = [df.velds.u.sel(time=ti)['time'].values for ti in Tlist]
        #Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list, T_list, df_extracted 
        #clean the memory
        gc.collect()  # Force garbage collection
    elif(param=="turbidity" or param=="avg_BS"):
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg  = np.mean(df.amp, axis = 0) 
        #Loop over the list of the time
        ########################
        P_list       = [np.nanmean(df_beam_avg.sel(time=ti), axis=0) for ti in Tlist]
        T_list       = [df_beam_avg.sel(time=ti)['time'].values for ti in Tlist]
        #Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list , T_list ,df_extracted 
        #clean the memory
        gc.collect()  # Force garbage collection
        ###########################
    elif(param=="dlnEcho"):
        #get the range form the data
        r           = df.range.data
        #make the average on all the 4 beams, NB the key word here is amp
        df_beam_avg = np.mean(df.amp, axis = 0)
#        Matrix2D    = np.asarray([df_beam_avg.sel(time=ti, range = r)  for ti in Tlist])
        Matrix2D   = xr.concat([df_beam_avg.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ####################################################
        T          = xr.concat([df_beam_avg.sel(time=ti)['time'] for ti in Tlist], dim="time")
#        #####################@
        DList       = np.copy(Matrix2D.values)
        #Free the memory for the large matrix
        del Matrix2D, df_beam_avg
        #clean the memory
        gc.collect()  # Force garbage collection
        #DList      = np.zeros(Matrix2D.values.shape, dtype = Matrix2D.values.dtype)
        for i in range(len(DList)):
            mean_val = np.mean(DList[i])
            DList[i] = [(item -mean_val)/mean_val for item in DList[i]]
        #DList = np.asarray(DList)
        if(average):
            ##%%%%%%%%%%%%% ##############
            #Return the Time and 1D perturbation in % (i.e. 100)
            return (T.values, np.mean(np.asarray(DList), axis=0) * 100) 
        ##############################
        #Return the Time and 2D perturbation in % (i.e. 100)
        return (T.values, np.asarray(DList) * 100)
        #################################################### 
    elif(param=="velocity2D_UP" or param=="velocity2D_DOWN"):
        #get the range form the data
        r         = df.range.data
        #get the Horizontal velocity
        U         = df.velds.U_mag
        #extract the 2D velocity
        Matrix2D  = xr.concat([U.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ###################################################
        T         = xr.concat([U.sel(time=ti)['time'] for ti in Tlist], dim="time")
        #####################@
        #Free the memory for the large matrix
        del U, r
        #clean the memory
        gc.collect()  # Force garbage collection
        #return the values
        return (T.values, Matrix2D.values)
    ##%%%%%%%% Extract the current direction %%%%%%%%%%%%%%%%%%
    elif(param=="veldir2D_UP" or  param=="veldir2D_DOWN"): 
        #get the range form the data
        r        = df.range.data
        #Get the velocity Direction in 2D
        U        = df.velds.U_dir
        #extract the 2D velocity
        Matrix2D = xr.concat([U.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ####################################################
        T        = xr.concat([U.sel(time=ti)['time'] for ti in Tlist], dim="time")
        #####################@
        #Free the memory for the large matrix
        del U, r
        #clean the memory
        gc.collect()  # Force garbage collection
        #return the values
        return (T.values, Matrix2D.values)
    ##%%%%%%%% Extract the current direction %%%%%%%%%%%%%%%%%%
    elif(param=="veldir1D_UP" or  param=="veldir1D_DOWN"): 
        #get the range form the data
        r        = df.range.data
        #Get the velocity Direction in 2D
        U        = df.velds.U_dir
        #extract the 2D velocity
        Matrix2D = xr.concat([U.sel(time= pd.to_datetime(ti).strftime("%Y-%m-%d"), range=r) for ti in sorted(Tlist)], dim="time")
        #T           = np.asarray([df_beam_avg.sel(time=ti)['time'].values for ti in Tlist])
        ####################################################
        T        = xr.concat([U.sel(time=ti)['time'] for ti in Tlist], dim="time")
        #####################@
        veldir1D = Matrix2D.mean(dim="range")
        #Free the memory for the large matrix
        del U, r
        #clean the memory
        gc.collect()  # Force garbage collection
        #return the values
        return (T.values, veldir1D)
    ################################################
    else:
        #Loop over the  list of the time an extract the desire parameter
        P_list       = [df[param].sel(time=ti) for ti in Tlist]
        T_list       = [df[param].sel(time=ti)['time'].values for ti in Tlist]
        # Flatten lists and create DataFrame
        df_extracted = pd.DataFrame({'time': np.hstack(T_list), 'value': np.hstack(P_list)})
        # Averaging duplicate timestamps
        df_extracted = df_extracted.groupby('time', as_index=False).mean()  # Averaging duplicate timestamps
        #################
        ####################################
        P            = df_extracted['value']
        T            = df_extracted['time'] 
        #################################
        del P_list, T_list, df_extracted 
    #Return the  value corresponding to the date of your choice
    return (T, P)

#########################
def Extract_database(df,  Tlist, param):
    #df is a dataframe and Date should be in the following format: YYYY-MM-DD
    df['date'] = df["DateTime"].dt.date
    dates_to_extract = [pd.to_datetime(year_month_day).date() for year_month_day in Tlist]
    #filtered for the selected days
    filtered = df[df['date'].isin(dates_to_extract)]
    #get the data for the selected days and the
    # Step 4: Extract wind_speed
    Data_Select = filtered[['DateTime', param]]
    #get the time and the data values
    Time_ALL    = Data_Select['DateTime'].to_numpy()
    Data_ALL    = Data_Select[param].to_numpy()
    #Free the memory
    del Data_Select, filtered
    #Return the  value corresponding to the date of your choice
    return (Time_ALL, Data_ALL)

######Add a new line in the dataset ###
def Add_Column(ds):
    ds['DateTime'] = pd.to_datetime(ds[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
    #Return the dataframe
    return ds




#Open the configuration file
with open("confdrose.yaml") as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)


#Grab all the parameters to plot
#Grab the Dictionary
PARAMSDICT          = Fig_params['PARAMSDICT']
#get the Meteo File
FILE_METEO          = Fig_params['FILE_METEO']
#get the ADCP File
ADCP_FILE_NAME_UP   = Fig_params['ADCP_FILE_NAME_UP']
ADCP_FILE_NAME_DOWN = Fig_params['ADCP_FILE_NAME_DOWN']
#wind rose diagram
wind_rose           = Fig_params['wind_rose']
#current rose diagram
#UP looking current
current_rose_up = Fig_params['current_rose_up']
#Down looking current
current_rose_down = Fig_params['current_rose_down']
#get the dataDate to plot
DictDate            = Fig_params['DictDate']
###################################
#transform into the list
DictDate            = list(DictDate)
#####################################

#read the Meteo database
d_mteo = pd.read_csv(FILE_METEO, delimiter ='\t')
#Add the of DateTime
d_mteo = Add_Column(d_mteo)


####### Read the ADCP FILES##############
try:
    df_down = dlfn.read(ADCP_FILE_NAME_DOWN)
    df_up   = dlfn.read(ADCP_FILE_NAME_UP)
except:
    print("print provide the ADCP files:  %s and %s  "%(ADCP_FILE_NAME_UP, ADCP_FILE_NAME_DOWN))
#plot_discharge(ax1,hours_of_day, data_day, starttime, Remove_frame)



#Grab only the parameter that will be plotted by arranging the PARAMSDICT and re-asigned a new avalues to PARAMSDICT
#PARAMSDICT = {k : PARAMSDICT[k]  for k in PARAMSDICT if(len(PARAMSDICT[k]) == 3)}
#print(PARAMSDICT)
#exit()
################################
#Set the number of subfigures to Zero
nfigs = 0
ii = 0
#Loop over the dictionary parameters
for i, key_param in zip(range(len(PARAMSDICT)), PARAMSDICT):
    #Check if the parameter is plotable
    if(Fig_params[key_param]):
        nfigs +=1
        #Add the index "i" in the list for numering the axis
        PARAMSDICT[key_param].append(ii)
        ii +=1
#change the values on the number of figure to plot by re-signing the len of PARAMSDICT
#nfigs  = len(PARAMSDICT)
#Grab the space between the figure
fig_space = Fig_params['fig_space']
#Grab the figure size
fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])
#Create the figure
#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size,constrained_layout=True)
#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
fig, axs  = plt.subplots(1,  nfigs, sharex = False, figsize = fig_size)

cmap = plt.get_cmap('viridis')  # This returns a colormap object

grey_shades = ['#d9d9d9', '#bdbdbd', '#969696', '#636363']
#bins        = np.linspace(0, 6, 4) 
bins        = np.linspace(0, 10, 4) 

if(wind_rose):
    #get the Dischage data for the corresponding time
    param         = "wind_rose"
    #get the wind speed data for the corresponding time
    time, data_spd    = Extract_database(d_mteo,  DictDate, "wind_speed")
    #extract the direction
    time, data_dir= Extract_database(d_mteo,  DictDate, "wind_direction")
    #get the parameter, the index of the corresponding axis, and the color
    title_xi , color, ix   = PARAMSDICT[param]
    #Plot the figure by calling the plotting function, plot_twinx
    #Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT) 
    #Remove placeholder axes so we can replace with WindroseAxes
    try:
        pos_ix = axs[ix].get_position()
        fig.delaxes(axs[ix])
    except:
        pos_ix = axs.get_position()
        fig.delaxes(axs)
    #Windrose for current
    ax_w   = WindroseAxes(fig, pos_ix)
    ######################
    # Set background to soft grey
    #ax_w.set_facecolor('#f0f0f0') 
    #Set grid appearance
    #ax_w.grid(True, color='white', linestyle='-', linewidth=1)
    # Define grey shades (lighter to darker)
    ########################
    fig.add_axes(ax_w)
    #ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='white')
    #ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='white',  antialiased=True, bins=bins,colors=grey_shades,  cmap=cmap)
    #ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='white',  antialiased=True, bins=bins,  cmap=cmap)
    ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='grey',  antialiased=True,  cmap=cmap)
    ax_w.set_title(title_xi)
    ax_w.set_legend(title="Speed (m/s)")

##Doubleched if we need to plot the Velocity
if(current_rose_up):
   #set the parameter
   param           = "current_rose_up"
   #get the horizontal velocity data for the corresponding time
   time, data_spd  = Extract_ADCPData(df_up , DictDate, "velocity_up")
   #####################
   time, data_dir  = Extract_ADCPData(df_up , DictDate, "veldir1D_UP")
   #get the parameter, the index of the corresponding axis, and the color
   title_xi , color, ix   = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #print(data_dir)
   #print(data_spd)
   #exit()
   #Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
   #Remove placeholder axes so we can replace with WindroseAxes
   try:
        pos_ix = axs[ix].get_position()
        fig.delaxes(axs[ix])
   except:
        pos_ix = axs.get_position()
        fig.delaxes(axs)
   #Windrose for current
   ax_curr = WindroseAxes(fig, pos_ix)
   fig.add_axes(ax_curr)
   #ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='white',  antialiased=True, bins=bins,  cmap=cmap)
   #ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='grey',  antialiased=True,  cmap=cmap)
   ax_w.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='grey', cmap=cmap)
   ax_curr.set_title(title_xi)
   ax_curr.set_legend(title="Speed (m/s)")
   #Plot the figure by calling the plotting function

if(current_rose_down):
   #set the parameter
   param           = "current_rose_down"
   #get the horizontal velocity data for the corresponding time
   time, data_spd  = Extract_ADCPData(df_down , DictDate, "velocity_down")
   #####################
   time, data_dir  = Extract_ADCPData(df_down , DictDate, "veldir1D_DOWN")
   #get the parameter, the index of the corresponding axis, and the color
   title_xi , color, ix   = PARAMSDICT[param]
   #Plot the figure by calling the plotting function
   #Plot_fig(axs[ix], time, data, DictDate, param, **PARAMSDICT)
   #Remove placeholder axes so we can replace with WindroseAxes
   try:
        pos_ix = axs[ix].get_position()
        fig.delaxes(axs[ix])
   except:
        pos_ix = axs.get_position()
        fig.delaxes(axs)
   #Windrose for current
   ax_curr = WindroseAxes(fig, pos_ix)
   fig.add_axes(ax_curr)
   #Define bins
   #rounded_max = np.ceil(np.max(data_spd))
   #bin_width = 2  # e.g. bins 0–2, 2–4, etc.
   #bins = np.arange(0, rounded_max + bin_width, bin_width)
   #n_bins = len(bins)   # Number of intervals
   #grey_map = cm.get_cmap('Greys', n_bins)
   #colors = [mcolors.to_hex(grey_map(i)) for i in range(n_bins)]
   #ax_curr.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='white')
   ax_curr.bar(data_dir, data_spd, normed=True, opening=0.8, edgecolor='grey', cmap=cmap)
   ax_curr.set_title(title_xi)
   ax_curr.set_legend(title="Speed (m/s)")

figname = "FigDirection-%s.pdf"%(list(DictDate)[0])

##############################################
#plt.xticks(fontsize= 18)
# Loop over each axis and set the tick font size
try:
    for ax in axs.flat:  # .flat flattens 2D array into 1D
        ax.tick_params(axis='x', labelsize=14)
except:
    pass
########################################
# Bring the window to the front (works with Qt5Agg or TkAgg backends)
#plt.subplots_adjust(hspace = 0.08)
#plt.subplots_adjust(hspace = fig_space)
#Align ylabel of the figure
fig.align_ylabels(axs)
#Save the figure
#fig.tight_layout()
#fig.canvas.manager.window.raise_()
#Save the figure
fig.savefig(figname, bbox_inches = 'tight')

