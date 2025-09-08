#########


from obspy import read_inventory
import matplotlib.pyplot as plt
import numpy as np

# Load the StationXML file

##########################################
Response_File= 'RESPOSE-FILE.XML/xj-fdsn.xml'
################################

# Load StationXML
inv = read_inventory(Response_File)
network = inv[0]
station = network[0]
channel = station[0]
response = channel.response

# Define parameters
sampling_rate = 250.0  # Update based on your instrument
nfft = 2 ** 14
frequencies = np.logspace(-2, 2, 1000)  # 0.01 to 100 Hz
t_samp = 1.0 / sampling_rate

# get_evalresp_response expects t_samp, nfft, and returns full FFT spectrum
resp, phase = response.get_evalresp_response(
    t_samp=t_samp,
    nfft=nfft,
    output="ACC"
)

# Determine the actual frequency axis
freq_axis = np.fft.rfftfreq(nfft, d=t_samp)

# Interpolate response at desired frequencies
response_interp = np.interp(frequencies, freq_axis, np.abs(resp))

# Convert to periods
periods = 1.0 / frequencies

# Plot response vs. period
plt.figure(figsize=(10, 6))
plt.semilogx(periods, response_interp)
plt.xlabel("Period [s]")
plt.ylabel("Amplitude [m/s]")
plt.title(f"Instrument Response for {channel.code} (vs. Period)")
plt.grid(True, which="both", ls="--")
plt.gca().invert_xaxis()
plt.tight_layout()

# Save to PDF
plt.savefig("response_vs_period.pdf", format="pdf", dpi=300)


