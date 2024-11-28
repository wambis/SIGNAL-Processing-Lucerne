#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
from glob import glob
#from math import log10, floor, ceil
import numpy as np
#from tabulate import tabulate
#for interpolation
#######################################
#######################
#import xlrd
#################################################
##############################################
#import matplotlib.dates as mdates
#######################################
#from matplotlib.colors import Normalize
############################
#import datetime as dtt 
#import xarray as xr
#import pandas as pd
#from mhkit import dolfyn
from mhkit import dolfyn as dlfn
#################
#import matplotlib.dates as mdates
#from mhkit.dolfyn.adp import api
#from pandas.core.common import flatten
##################################




import pandas as pd
from dolfyn import read
import datetime


import xarray as xr

#print(dir(xr.Dataset))
#
#print(xr.__version__)

#print(pd.__version__)
#exit()



class ADCPDataExtractor:
   def __init__(self, FilePath):
        self.FilePath = FilePath
        self.data     = None

   def load_data(self):
       """ Load ADCP data from the FilePath """
       self.data = read(self.FilePath) 

   def Average(self, n_bin=1):
    #To average the data into time bins, we need the parmeter
    #Define functio
    #1) n_bin which represents the number of data points in each ensemble,
    #2) here our enseblme happen every 600s and  the frequency is fs = 0.001668, so n_bin = fs* 600 = 1

    #n_bin = 1 is the value by default
    #Set the average tool
    avg_tool = dlfn.VelBinner(n_bin = n_bin, fs = self.data.fs)
    #Perform the averge
    self.data   = avg_tool.bin_average(self.data)
    #Compute the average of the horizontal velocity known as U_mag and it direction known as U_dir
    #return self.data


   def extract_one_day_data(self, target_date, param, nrange= None, save_in_file=None):
        #if(nrange== None):
        #    r        = df.range.data
        #else:
        #    #r        = df.range.data[:nrange]
        #    r        = self.data.range[:nrange]
        ############################################
        if not isinstance(target_date, str):
            print("print the time should be a string on the format: YYYY-MM-DD")
            exit(1)
        #Check the paramters
        if(param=='U_mag'):
                #get the vertical velocity of  Horizontal velocity water
                try:
                    data     = np.nanmean(self.data.velds.U_mag.sel(time = target_date), axis =0)
                    time     = self.data.velds.U_mag.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='u'):
                #get the vertical velocity of  Northward water
                try:
                    data = np.nanmean(self.data.velds.u.sel(time = target_date), axis =0)
                    time = self.data.velds.u.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='v'):
                #get the vertical velocity of  Easternward water
                try:
                    data  = np.nanmean(self.data.velds.v.sel(time = target_date), axis =0)
                    time  = self.data.velds.v.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='w'):
                #get the vertical velocity of  the water
                try:
                    data  = np.nanmean(self.data.velds.w.sel(time = target_date), axis =0)
                    time = self.data.velds.w.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 

        elif(param=='U_dir'):
                #get  the water Direction measureb by the ADCP
                try:
                    data = np.nanmean(self.data.velds.U_dir.sel(time = target_date), axis =0)
                    time = self.data.velds.U_dir.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 

        elif(param=="amp"):
                # get the backscatter from the water
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(self.data.amp, axis = 0) 
                try:
                    data  = np.nanmean(df_beam_avg.sel(time = target_date), axis =0)
                    time  = df_beam_avg.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
            #Loop over the list of the time
        else:
                try:
                    data  = self.data[param].sel(time = target_date)
                    data  = data.data
                    time  = self.data[param].sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        #Check if we need to save into a file
        if(save_in_file):
            #Create a DataFrame
            dff = pd.DataFrame({'time': time, param: data})
            #save into a file
            dff.to_csv("DataExtracted.csv", index=False)

        return (data, time)


   def extract_one_day_data2D(self, target_date, param, nrange= None):
        if(nrange== None):
            r        = df.range.data
        else:
            #r        = df.range.data[:nrange]
            r        = self.data.range[:nrange]
        ############################################
        if not isinstance(target_date, str):
            print("print the time should be a string on the format: YYYY-MM-DD")
            exit(1)
        #Check the paramters
        if(param=='U_mag'):
                #get the vertical velocity of  Horizontal velocity water
                try:
                    Data     = self.data.velds.U_mag.sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data.velds.U_mag.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='u'):
                #get the vertical velocity of  Northward water
                try:
                    Data     = self.data.velds.u.sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data.velds.u.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='v'):
                #get the vertical velocity of  Easternward water
                try:
                    Data     = self.data.velds.v.sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data.velds.v.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        elif(param=='w'):
                #get the vertical velocity of  the water
                try:
                    Data     = self.data.velds.w.sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data.velds.w.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 

        elif(param=='U_dir'):
                #get  the water Direction measureb by the ADCP
                try:
                    Data     = self.data.velds.U_dir.sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data.velds.U_dir.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 

        elif(param=="amp"):
                #get the backscatter from the water
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(self.data.amp, axis = 0) 
                try:
                    Data  = df_beam_avg.sel(time = target_date, range = r)
                    data  = Data.data
                    time  = df_beam_avg.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
            #Loop over the list of the time
        else:
                try:
                    Data     = self.data[param].sel(time = target_date, range = r)
                    data     = Data.data
                    time     = self.data[param].sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
        ##Check if we need to save into a file
        #if(save_in_file):
        #    #Create a DataFrame
        #    dff = pd.DataFrame(data,  columns=[param, param])
        #    #Add the time as the first column
        #    dff.insert(0, "time", time.T)              
        #    #save into a file
        #    dff.to_csv("DataExtracted.csv", index=False)

        return (data, time)


   def extract_several_days_data(self, target_date, param, nrange= None, save_in_file = None):
        ############################################
        if not isinstance(target_date, list):
            print("print the time should be a List of on the format date in the format: YYYY-MM-DD")
            exit(1)
        if(param=='U_mag'):
            #get the vertical velocity of  Horizontal velocity water
            try:
                data  = np.concatenate([np.nanmean(self.data.velds.U_mag.sel(time = ti), axis =0) for ti in target_date])
                time  = np.concatenate([self.data.velds.U_mag.sel(time = ti)['time'] for ti in target_date])
            except:
                raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 
        elif(param=='u'):
                #get the vertical velocity of  Northward water
            try:
                data  = np.concatenate([np.nanmean(self.data.velds.u.sel(time = ti), axis =0) for ti in target_date])
                time  = np.concatenate([self.data.velds.u.sel(time = ti)['time'] for ti in target_date])
            except:
                raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 
        elif(param=='v'):
                #get the vertical velocity of  Easternward water
            try:
                data  = np.concatenate([np.nanmean(self.data.velds.v.sel(time = ti), axis =0) for ti in target_date])
                time  = np.concatenate([self.data.velds.v.sel(time = ti)['time'] for ti in target_date])
            except:
                raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 
        elif(param=='w'):
                #get the vertical velocity of  the water
            try:
                data  = np.concatenate([np.nanmean(self.data.velds.w.sel(time = ti), axis =0) for ti in target_date])
                time  = np.concatenate([self.data.velds.w.sel(time = ti)['time'] for ti in target_date])
            except:
                raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 

        elif(param=='U_dir'):
                #get  the water Direction measureb by the ADCP
            try:
                data  = np.concatenate([np.nanmean(self.data.velds.U_dir.sel(time = ti), axis =0) for ti in target_date])
                time  = np.concatenate([self.data.velds.U_dir.sel(time = ti)['time'] for ti in target_date])
            except:
                raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 

        elif(param=="amp"):
                # get the backscatter from the water
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(self.data.amp, axis = 0) 
                try:
                      #Loop over the list of the time
                      data  = np.concatenate([np.nanmean(df_beam_avg.sel(time = ti), axis =0) for ti in target_date])
                      time  = np.concatenate([df_beam_avg.sel(time = ti)['time'] for ti in target_date])
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 
        else:
                try:
                    #Loop over the  list of the time an extract the desire parameter
                    data  = np.concatenate([ self.data[param].sel(time = ti) for ti in target_date])
                    time  = np.concatenate([self.data[param].sel(time = ti)['time'] for ti in target_date])
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, ' , '.join(target_date))) 
        #Check if we need to save into a file
        if(save_in_file):
            #Create a DataFrame
            dff = pd.DataFrame({'time': time, param: data})
            #save into a file
            dff.to_csv("DataExtracted.csv", index=False)

        return (data, time)

    #Check if the desired velocity is the Horizontal velocity

   def extract_day_data(self, target_date):
      """ 
      Extract the data for a specific day.
       
      Parameter:
      target_date (str) : The target date is in 'YYYY-MM-DD' format

      Return:
      pd.DataFrame containing the ADCP Data for a specific day.
      """
      if self.data is None:
           raise ValueError('Data has not been load, please call load_data() first.') 

      #Convert the target date to datetime object

      target_date = pd.to_datetime(target_date)
      #pd.datetime Converts date strings into pandas datetime objects

      #Filter the data for the target date
      start_time = target_date
      end_time   = target_date + pd.Timedelta(days=1)
      #Assuming the data has a `time` column
      day_data  = self.data[ (self.data['time'] >= start_time) & (self.data['time'] < end_time) ]

      #return the extracted data
      return day_data


#Usage  example

if __name__ == "__main__":
    File_path = "DATASELECT/MUO2_UP_11_23.000"

    target_date = '2023-10-23'  # Example target date
    #target_date = ['2023-10-22', '2023-10-23']  # Example target date
    #get the extractor object
    extractor = ADCPDataExtractor(File_path)
    #############
    extractor.load_data()
    #print(extractor.Average(n_bin=1))
    extractor.Average(n_bin=1)
    #P, T = extractor.extract_data(target_date, nrange=15)
    P, T = extractor.extract_one_day_data(target_date, "pressure", nrange=15, save_in_file= True)
    #P, T  = extractor.extract_several_days_data(target_date, "pressure", save_in_file= True)
    #P, T  = extractor.extract_several_days_data(target_date, "prcnt_gd", save_in_file= True)
    #P, T  = extractor.extract_several_days_data(target_date, "U_mag", nrange=15,save_in_file= True)
    #P, T  = extractor.extract_one_day_data2D(target_date, "pressure", nrange=15) 
    #print(day_data)
    print(P)
    print(P.shape)
    print(T.shape)

