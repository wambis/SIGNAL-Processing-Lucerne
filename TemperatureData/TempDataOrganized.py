#!/usr/bin/env python
#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
from glob import glob
import pandas as pd



Dirs = glob("AUTO_*") 

#print(Dirs)
#exit()

filecode= "conf-Thermistors-2024-06-07.yaml"


with open(filecode) as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    FILE_params = yaml.load(Fym, Loader=SafeLoader)

#print(FILE_params)
#exit()
#Open a dictionary for data
DataDict =dict()

for FILE in Dirs:
    if os.path.isdir(FILE):
        #print(f"'{FILE}' is a directory!")
        fcode = FILE.split("_")[1]
        #ftemp = glob("%s/AUTO_200566_20240607_0941_data.txt
        try:
            ftemp = glob("%s/AUTO_%s*_data.txt"%(FILE,str(fcode)))
        except:
            print("The file corresponding to the code %s doest not exits!"%(str(fcode))) 
            exit()
        if(len(ftemp) > 1):
            print("The code %s corresponds to two files %s"%(fcode, ftemp))
            exit()
        else:
            ftemp = ftemp[0]
        #####
        try:
            depth= FILE_params["code_depth"][fcode]
        except:
            print("this file code %s does not correspond to any depth"%(str(fcode)))
            exit()

        #print(fcode)
        #Read CSV into DataFrame
        #df = pd.read_csv(ftemp, parse_dates=['Time'])
        df = pd.read_csv(ftemp)
        #Convert to dictionary
        #data_dict = dict(zip(df['Time'].astype(str), df['Temperature']))
        data_dict       = dict(zip(df['Time'], df['Temperature']))
        #save in the the Data Dictionnary
        DataDict[depth] = data_dict
        #print(DataDict)
        #exit()

## Specify the file name to save the temperature data
file_name = "DataTmeperature.json"

# Write data to a JSON file
with open(file_name, "w") as json_file:
    json.dump(DataDict, json_file, indent=4)

print(f"JSON file '{file_name}' created successfully!")
