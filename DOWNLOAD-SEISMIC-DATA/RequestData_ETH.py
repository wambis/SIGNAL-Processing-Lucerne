#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from datetime import datetime, timedelta
# Initialize the client (IRIS or another FDSN service supporting CH network)
#client = Client("IRIS")  # Replace "IRIS" with the appropriate provider if needed
client = Client("ETH")  # Replace "ETH" with the appropriate provider if needed
# Define parameters
network = "CH"
station = "MUO"
location = "*"  # Use '*' to include all locations or specify if known
#channel = "*"   # Use '*' to include all channels or specify if known
channel = "BH*"   # Use '*' Broad band High gain channel
#channel = "HH*"   # High sampling High gain broad band channel

#######################
#start_day = "2023-10-07T00:00:00"
start_day = "2023-10-07"
#set the number of the days your wish to download the data
number_of_days = 1

#start date
start_date = datetime.strptime("%sT00:00:00"%(start_day), "%Y-%m-%dT%H:%M:%S")
#end_date
end_date   = start_date + timedelta(days = number_of_days) 
#get the starttime
start_time = UTCDateTime(start_date)  # Start of the day
#get the endtime
end_time   = UTCDateTime(end_date)    # End of the day
#get the last day
#year_month_day = "%s-%s-%s" % (end_date.year, end_date.month, end_date.day) 
Iyear_month_day = "%s-%s-%s" % (end_time.year, end_time.month, end_time.day) 
Lyear_month_day = "%s-%s-%s" % (start_time.year, start_time.month, start_time.day) 

#if(Iyear_month_day != Lyear_month_day):
if(number_of_days > 1):
    #set the output MiniSEED file name
    mseed_filename = "MUO.CH_%s_to_%s.mseed"%(start_day, Lyear_month_day)
    #set the output of the xml file name
    xml_filename   = "MUO.CH_%s_to_%s.xml"%(start_day, Lyear_month_day)
else:
    #set the output MiniSEED file name
    mseed_filename = "MUO.CH_%s.mseed"%(start_day)
    #set the output of the xml file name
    xml_filename   = "MUO.CH_%s.xml"%(start_day)
# Download MiniSEED data
print("Downloading MiniSEED data...")
st = client.get_waveforms(network, station, location, channel, start_time, end_time)
#mseed_filename = "MUO_2024-01-01.mseed"
#Write into the MiniSEED file
st.write(mseed_filename, format="MSEED")
print(f"MiniSEED data saved to {mseed_filename}")
#Download StationXML metadata
print("Downloading StationXML metadata...")
inventory    = client.get_stations(network=network, station=station, level="response")
#xml_filename = "MUO_2024-01-01.xml"
inventory.write(xml_filename, format="STATIONXML")
print(f"StationXML metadata saved to {xml_filename}")
