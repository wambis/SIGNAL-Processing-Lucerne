#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
from glob import glob
#from math import log10, floor, ceil
import numpy as np
from mhkit import dolfyn as dlfn
import matplotlib.dates as mdates
import matplotlib.pyplot as plt 
from obspy import read, read_inventory
from obspy.core import UTCDateTime
############################
import pandas as pd
import datetime
import xarray as xr


""" 
Class writen by @Mathurin Wamba@      
University of Bern, December, 2024
This packages is used for ADCP data extraction as well as
Meteorological extraction and the input files haves to be 
Organized as those provide as example.
The user has to provide the 
Parameter: U_mag, u, v, w, amp, U_dir, pressure etc...
target_date (str) : The target date is in 'YYYY-MM-DD' format
and the package will return
the Data containing the ADCP Data for a specific or several days.
"""


class ADCPDataExtractor:
   def __init__(self, FilePath):
        self.FilePath = FilePath
        self.data     = None

   def load_data(self):
       """ Load ADCP data from the FilePath """
       self.data = dlfn.read(self.FilePath) 

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


   def extract_one_day_data(self, target_date, param, nbin= None, starttime= None, endtime=None, save_in_file=None):
       #grab the range
        if(nbin== None):
            r        = self.data.range
        else:
            try:
                if(nbin >= 1):
                    r   = self.data.range[nbin -1 : nbin]
                else:
                    r   = self.data.range[:nbin]
            except:
                raise ValueError("The bin number %i seems not to exist, please check the bin number "%(nbin)) 
        ############################################
        if not isinstance(target_date, str):
            print("print the time should be a string on the format: YYYY-MM-DD")
            exit(1)
        #Check the paramters
        if(param=='U_mag'):
                #get the vertical velocity of  Horizontal velocity water
                try:
                    data      = np.nanmean(self.data.velds.U_mag.sel(time = target_date, range= r), axis =0)
                    time      = self.data.velds.U_mag.sel(time = target_date)['time'] 
                    print(data)
                    DataFrame = pd.DataFrame(data, index=time)
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #print(data)
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        elif(param=='u'):
                #get the vertical velocity of  Northward water
                try:
                    data = np.nanmean(self.data.velds.u.sel(time = target_date, range= r), axis =0)
                    time = self.data.velds.u.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        elif(param=='v'):
                #get the vertical velocity of  Easternward water
                try:
                    data  = np.nanmean(self.data.velds.v.sel(time = target_date, range= r), axis =0)
                    time  = self.data.velds.v.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        elif(param=='w'):
                #get the vertical velocity of  the water
                try:
                    data  = np.nanmean(self.data.velds.w.sel(time = target_date, range= r), axis =0)
                    time = self.data.velds.w.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        elif(param=='U_dir'):
                #get  the water Direction measureb by the ADCP
                try:
                    data = np.nanmean(self.data.velds.U_dir.sel(time = target_date, range= r), axis =0)
                    time = self.data.velds.U_dir.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        elif(param=="amp"):
                # get the backscatter from the water
                #make the average on all the 4 beams, NB the key word here is amp
                df_beam_avg = np.mean(self.data.amp, axis = 0) 
                try:
                    data  = np.nanmean(df_beam_avg.sel(time = target_date, range= r), axis =0)
                    time  = df_beam_avg.sel(time = target_date)['time'] 
                except:
                    raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        #Loop over the list of the time
        else:
                try:
                    data  = self.data[param].sel(time = target_date)
                    data  = data.data
                    time  = self.data[param].sel(time = target_date)['time'] 
                except:
                    try:
                        data = np.nanmean(self.data[param].sel(time = target_date, range= r), axis =0)
                        time = self.data[param].sel(time = target_date)['time'] 
                    except:
                        raise ValueError("The parameter %s or the date %s seems not to exist, please check them"%(param, target_date)) 
                #check the extraction need to be at the specific periode of a day
                if(starttime !=None and endtime !=None):
                    t1       = '%s %s'%(target_date, starttime)
                    t2       = '%s %s'%(target_date, endtime)
                    try:
                        df_intrv = DataFrame[t1 : t2]
                    except:
                        raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, target_date)) 
                    #Convert DatetimeIndex to matplotlib dates 
                    time     = mdates.date2num(df_intrv.index)
                    #get the data
                    data     = df_intrv.iloc[:, 0]
                    #Convert pandas Series to NumPy array 
                    data     = data.to_numpy()
        #Check if we need to save into a file
        if(save_in_file):
            #Create a DataFrame
            dff = pd.DataFrame({'time': time, param: data})
            #save into a file
            dff.to_csv("DataExtracted.csv", index=False)

        return (time, data)


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

        return (time, data)


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

        return (time, data)

#A class to extract the meteorological Data
class MeteoDataExtractor:
   def __init__(self, FilePath):
        self.FilePath = FilePath
        self.data     = None

   #read the Meteo database
   def load_data(self):
       """ Load meteorological data from the FilePath """
       self.data = pd.read_csv(self.FilePath, delimiter ='\t') 
       #Add the of DateTime

   def AddDateTimeColumn(self, dparam, save_in_file = None):
        """ 
        Add DateTime as first column in the dataset and 
        clean the data itself 
        The Data format should be in the form Year (YR), Month (MO), Day (DA), Hours (HH)
        Minutes (MM), Second (SS) 
        The Unnamed column will be removed from the Dataset. 
        """
        #self.data['DateTime'] = pd.to_datetime(self.data[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
        self.data[dparam] = pd.to_datetime(self.data[['YR', 'MO', 'DA', 'HH', 'MM', 'SS']].astype(str).agg(' '.join, axis=1), format='%y %m %d %H %M %S')
        ##########################
        #Drop the original four columns 
        self.data        = self.data.drop(columns=['YR', 'MO', 'DA', 'HH', 'MM', 'SS'])
        #Drop columns with names that contain 'Unnamed' 
        self.data        = self.data.loc[:, ~self.data.columns.str.contains('^Unnamed')]
        #Reorder the columns to make sure 'col1' is the first column 
        cols             = [dparam] + [col for col in self.data.columns if col != dparam]
        self.data        = self.data[cols] 
        #Set the datetime column as the index 
        self.data.set_index(dparam,  inplace=True)
        if(save_in_file):
            #Save to CSV file 
            self.data.to_csv('data.csv', sep='\t',index=True)
        else:
            #Return the dataframe
            return self.data

   def extract_one_day_mdata(self, date_of_day, param, starttime= None, endtime=None, save_in_file=None):
       #Extract data for a specific day in the format YYYY-MM-DD
       day_data       = self.data.loc[date_of_day]
       day_data_param = day_data[param].to_numpy()  
       #Convert DatetimeIndex to matplotlib dates 
       time           = mdates.date2num(day_data.index)
       #########################
       #check the extraction need to be at the specific periode of a day
       if(starttime!=None and endtime!=None):
           #set the time interval to be extracted
           t1       = '%s %s'%(date_of_day, starttime)
           t2       = '%s %s'%(date_of_day, endtime)
           try:
               #Extract the dta
               df_in = day_data[t1 : t2]
           except:
               raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, date_of_day)) 
           day_data_param = df_in[param].to_numpy() 
           time           = mdates.date2num(df_in.index)
       #Return the value
       return (time, day_data_param)

#A class to extract the meteorological Data
class SeismometerInfoExtractor:
   def __init__(self, FilePath):
        self.FilePath = FilePath
        self.data     = None

   #read the Meteo database
   def load_data(self):
       """ Load seismometer information from the FilePath """
       self.data = pd.read_csv(self.FilePath) 
       #Originals columns
       ocols     = [it for it in self.data.columns]
       #Arrange the data structure
       cols      = [it.split('[')[0].strip() for it in self.data.columns]
       cols      = ['%s-%s'%(it.split()[0], it.split()[1]) if 'Battery' in it else it  for it in cols]
       #change the colomn name by removig °, [ and switching Battery Voltage to Battery-Voltage
       #Rename columns
       try:
           self.data.rename(columns={ocols[1] : cols[1], ocols[2]: cols[2], ocols[3]: cols[3]}, inplace=True)
       except:
               raise ValueError("Not able to renames the columns %s  %s  %s"%(cols[1], cols[2], cols[3])) 
       # Remove "UTC" from the datetime strings 
       #self.data['Time'] = self.data['Time'].str.replace(' UTC', '')
       #return self.data

       #Add the of DateTime
   def extract_one_day_sinfo(self, date_of_day, param, starttime= None, endtime=None, save_in_file=None):
       try:
           #Remove "UTC" from the datetime strings 
           self.data['Time'] = self.data['Time'].str.replace(' UTC', '')
           #Convert the 'datetime' column to datetime format, so that we can extract the daily data 
           self.data['Time'] = pd.to_datetime(self.data['Time'])
       except:
           pass
       #Set the datetime column as the index 
       cols           = [it for it in self.data.columns]
       self.data.set_index(cols[0],  inplace=True)
       print(self.data)
       #exit()
       #Extract data for a specific day in the format YYYY-MM-DD
       data_day       = self.data.loc[date_of_day]
       data_day_param = data_day[param].to_numpy()  
       #Convert DatetimeIndex to matplotlib dates 
       time           = mdates.date2num(data_day.index)
       #########################
       #check the extraction need to be at the specific periode of a day
       if(starttime!=None and endtime!=None):
           #set the time interval to be extracted
           t1         = '%s %s'%(date_of_day, starttime)
           t2         = '%s %s'%(date_of_day, endtime)
           try:
               #Extract the dta
               df_in= data_day[t1 : t2]
           except:
               raise ValueError("Not able to extract the data between %s and  %s of the day: %s"%(t1, t2, date_of_day)) 
           data_day_param = df_in[param].to_numpy() 
           time           = mdates.date2num(df_in.index)
       #Return the value
       return (time, data_day_param)



def read_mseed(Date_of_Day, mseedFile, xmlFile, bandpass, Pressure=None, start_time=None, end_time=None):
    #get the stream by reading the mseed file
    st      = read(mseedFile)
    #read the inventory
    inv     = read_inventory(xmlFile)
    #prefiltering the signal
    #Very Good match with the discharge of the day with this event  XJ_MUS01_HH2-2018-12-24.mseed 
    dT      = '%s s'%(bandpass)
    #Periods List
    pds     = bandpass.split("-") 
    f1      = 1./float(pds[-1]); f2 = 1./float(pds[-2]) ; 
    ####################################################
    f3      = 1./float(pds[1]);  f4 = 1./float(pds[0])
    pre_filt = [f1, f2, f3, f4]
    if(Pressure):
        st.remove_response(inventory=inv, pre_filt=pre_filt, output="DEF",water_level=60)
    else:
        #Remove instrument response
        st.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60)
        #st.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL",water_level=60)
    #get the trace from the stram
    for tr in st:
        #get the parameters of the trace
        network       = tr.stats.network
        station       = tr.stats.station
        channel       = tr.stats.channel
        stime         = str(tr.stats.starttime)
        endtime       = str(tr.stats.endtime)
        #remove mean and trend on the trace
        tr.detrend(type='demean')
        Title         = "%s  %s"%(network, station)
        ##Resampling the data
        tr.resample(sampling_rate = 1.0, window = 'hann', no_filter = True)
        #################################
        if(start_time !=None and end_time != None):
            #Define the start and end times
            t1         = "%sT%s"%(Date_of_Day, start_time)
            t2         = "%sT%s"%(Date_of_Day, end_time)
            start_time = UTCDateTime(t1)
            end_time   = UTCDateTime(t2)
            dtt            = end_time - start_time
            #Trim the trace between the specified times
            tr         = tr.trim(starttime=start_time, endtime=end_time)
            #tr         = tr.trim(start_time, start_time + dtt)
            data       = tr.data
            date_rng   = pd.date_range(start= t1, end = t2, freq='1s')
           #sampling_rate = tr.stats.sampling_rate
        else:
            data       = tr.data
            time       = tr.times()
            #The sampling rate is 1 second
            date_rng  = pd.date_range(start= Date_of_Day, freq='1s', periods= data.size)
        #change the date_rng to numpy
        time      = date_rng.to_numpy()
    return(time, tr, Title)
#Usage  example

if __name__ == "__main__":
    File_path       = "DATASELECT/MUO2_UP_11_23.000"
    FileMeteopath   = "DATASELECT/DATA_ARRANGE_ALT.csv"
    Filesispath     = "DATASELECT/MUS03.aux"
    ###################################################
    mseedFile       = "STATION03/XJ_MUS03_HH2-2023-10-23.mseed" 
    xmlFile         = "RESPOSE-FILE.XML/xj-fdsn.xml"
    #############################################
    target_date     = '2023-10-23'  # Example target date
    bandpass        = '0.02-0.022-200-1000'
    #target_date = ['2023-10-22', '2023-10-23']  # Example target date
    #get the extractor object
#    extractor   = ADCPDataExtractor(File_path)
#    #############
#    extractor.load_data()
#    #print(extractor.Average(n_bin=1))
#    extractor.Average(n_bin=1)
#    #T, P = extractor.extract_one_day_data(target_date, "U_mag", starttime= '10:00', endtime='14:30', save_in_file= False)
#    T, P = extractor.extract_one_day_data(target_date, "U_mag", nbin= 10, save_in_file= False)
#    #T, P = extractor.extract_one_day_data(target_date, "pressure", nrange=15, save_in_file= True)
#    #T, P  = extractor.extract_several_days_data(target_date, "pressure", save_in_file= True)
#    #T, P  = extractor.extract_several_days_data(target_date, "prcnt_gd", save_in_file= True)
#    #T, P  = extractor.extract_several_days_data(target_date, "U_mag", nrange=15,save_in_file= True)
#    #T, P  = extractor.extract_one_day_data2D(target_date, "pressure", nrange=15) 
#    #print(day_data)
#    #print(T)
#    print(P.shape)
#    print(type(T))
#    #####################################################
    fig, ax    = plt.subplots(1, 1, sharex = False, figsize = (12,11))
#    #plot the frequency as function of time
#    #ax.plot(times, freqs, lw='2', c='r')
#    ax.plot(T, P, lw=1.0, c='r')
    #time, tr, Title =  read_mseed(target_date, mseedFile, xmlFile, bandpass, Pressure=None, start_time='21:00', end_time='23:30')
    #time, tr, Title =  read_mseed(target_date, mseedFile, xmlFile, bandpass, Pressure=None, start_time='17:30', end_time='21:00')
    time, tr, Title =  read_mseed(target_date, mseedFile, xmlFile, bandpass, Pressure=None, start_time='14:00', end_time='17:00')
    ax.plot(time, tr.data, lw=0.5, c='k', alpha = 0.5)
    print(time)
    print(type(time)) 
    #exit()

#    ###################################################
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    figname = "Freqs_TEST_PLOTS.png"
    plt.savefig(figname, bbox_inches = 'tight')
#    Mextractor = MeteoDataExtractor(FileMeteopath)
#    Mextractor.load_data() 
#    #Add the Data column and clean the data
#    Mextractor.AddDateTimeColumn('DateTime',save_in_file = True)
#    t, d = Mextractor.extract_one_day_mdata(target_date, 'wind_speed',starttime='10:00', endtime='14:30')
#    print(d)
#    print('*******' * 20) 
#    print(t)
#    print(len(t), len(d))
#    Seis  = SeismometerInfoExtractor(Filesispath)
#    Seis.load_data() 
#    #t, d = Seis.extract_one_day_sinfo(target_date, 'Temperature')
#    #t, d = Seis.extract_one_day_sinfo("2023-11-29", "Temperature")
#    t, d = Seis.extract_one_day_sinfo("2023-11-29", "Temperature", starttime='10:00', endtime='14:30')
#    print(d)

