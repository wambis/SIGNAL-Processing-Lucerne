#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
# Initialize the client (IRIS or another FDSN service supporting CH network)
#client = Client("IRIS")  # Replace "IRIS" with the appropriate provider if needed
client = Client("ETH")  # Replace "ETH" with the appropriate provider if needed
# Define parameters
network = "CH"
station = "MUO"
location = "*"  # Use '*' to include all locations or specify if known
#channel = "*"   # Use '*' to include all channels or specify if known
channel = "BH*"   # Use '*' to include all channels or specify if known
#channel = "HH*"   # Use '*' to include all channels or specify if known

#start_time = UTCDateTime("2023-10-23T00:00:00")  # Start of the day
#end_time = UTCDateTime("2023-10-24T00:00:00")    # End of the day

start_time = UTCDateTime("2023-10-07T00:00:00")  # Start of the day
end_time = UTCDateTime("2023-10-08T00:00:00")    # End of the day
# Download MiniSEED data
print("Downloading MiniSEED data...")
st = client.get_waveforms(network, station, location, channel, start_time, end_time)
#mseed_filename = "MUO_2023-10-23.mseed"
#mseed_filename = "MUO_2024-01-01.mseed"
mseed_filename = "MUO_2023-10-07.mseed"
st.write(mseed_filename, format="MSEED")
print(f"MiniSEED data saved to {mseed_filename}")

# Download StationXML metadata
print("Downloading StationXML metadata...")
inventory    = client.get_stations(network=network, station=station, level="response")
#xml_filename = "MUO_2023-10-23.xml"
#xml_filename = "MUO_2024-01-01.xml"
xml_filename = "MUO_2023-10-07.xml"
inventory.write(xml_filename, format="STATIONXML")
print(f"StationXML metadata saved to {xml_filename}")
