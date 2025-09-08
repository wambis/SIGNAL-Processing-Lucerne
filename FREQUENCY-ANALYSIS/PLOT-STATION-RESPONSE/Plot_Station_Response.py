##########
from obspy import read_inventory
import numpy as np
import matplotlib.pyplot as plt




##########################################
Response_File= 'RESPOSE-FILE.XML/xj-fdsn.xml'
# Load StationXML
inv = read_inventory(Response_File)

#Get response from first channel
network = inv[0]
station = network[0]
channel = station[0]
response = channel.response
#XJ.MUS05
print(network, station, channel)
#print(response)
#exit()
# Plot instrument response
#fig = response.plot(min_freq=0.01, output="VEL")
fig = response.plot(min_freq=0.001, output="VEL")

# Save to PDF
fig.savefig("station_response.pdf", format="pdf")

plt.show()

