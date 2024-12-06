#!/usr/bin/env python

#Plot the mseed seismic data files

import obspy

from obspy import read, read_inventory

import matplotlib.pyplot as plt

from glob import glob

import numpy as np

#from Ipython.display import display

import pandas as pd


def Change_time_notation(TimeFormat, year):
            #grab the time parameter
            if(TimeFormat.startswith(str(year))):
                tt = TimeFormat.replace(" ", ":")
            else:
                day_month_year, times = TimeFormat.split(" ")
                day, month, year      = day_month_year.split("-")
                tt = "%s-%s-%s:%s" % (year, month, day,times)
            return tt







#Name of the file to be read
fname       = "Discharge_2023.csv"

#load the data frame
df          = pd.read_csv(fname, skiprows=0, sep =";", encoding="us-ascii")

################################################################
NEW_TIMES = []

for k in df.keys():
    if(k=="DateTime"):
         for t in df["DateTime"]:
             if(t.startswith("2023")):
                tt = t.replace(" ", ":")
                NEW_TIMES.append(tt)
             else:
                day_month_year, times = t.split(" ")
                day, month, year      = day_month_year.split("-")
                tt = "%s-%s-%s:%s" % (year, month, day,times)
                NEW_TIMES.append(tt)

df['DateTime'] = NEW_TIMES

#Save the dataframe into a file
df.to_csv("Discharge_2023_Wamba.csv")




########################################
#stn_name    = df['Stationsname'] 
#stn_num     = df['Stationsnummer'] 
#params      = df['Parameter'] 
#########################################
#time_series = df['Zeitreihe']
#param_unit  = df['Parametereinheit'] 
#water_body  = df['Gewasser'] 
###################################
#DataTime    = df['DateTime'] 
#time_of_occ = df['Zeitpunkt_des_Auftretens'] 
#disch_value = df['Wert'] 
#relase_stat = df['Freigabestatus']

#Create a dictionary to keep all the data
#Dict        = {k : [] for k in df.keys()} 
####################

#for sname, snum, p, time_s, p_unit, wb, DT, tocc, dvalue, rst in zip(stn_name, stn_num, params,time_series,param_unit, water_body, DataTime, time_of_occ, disch_value, relase_stat):
#            #grab the time parameter
#            if(DT.startswith("2023")):
#                tt = DT.replace(" ", ":")
#            else:
#                day_month_year, times = DT.split(" ")
#                day, month, year      = day_month_year.split("-")
#                tt = "%s-%s-%s:%s" % (year, month, day,times)
