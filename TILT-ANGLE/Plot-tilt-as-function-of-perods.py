# ============================================================
# DayNoise / ATaCR Tilt Analysis
#
# Plot spectra as a function of PERIOD (s)
# with a TRUE logarithmic period axis.
#
# DayNoise calculations are performed in frequency (Hz).
# After DayNoise creates the figure, spectral curves are
# converted from frequency to period:
#
#                   T = 1 / f
#
# ============================================================


# ============================================================
# 1. Imports
# ============================================================

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from obspy import read
from obstools.atacr import DayNoise


# ============================================================
# 2. Define input/output paths
# ============================================================

base_path = Path(
    #"LucerneData/muota-SingSta-05/mseed/2024/01/13"
    "LucerneData/muota-SingSta-03/mseed/2023/09/22"
    #"LucerneData/muota-SingSta-05/mseed/2024/03/29"
)

output_dir = Path("./figures")

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 3. Define miniSEED files
# ============================================================
#file_h1 = base_path / "XJ_MUS05_HH2-2024-03-29.mseed"
#file_h2 = base_path / "XJ_MUS05_HH3-2024-03-29.mseed"
#file_z  = base_path / "XJ_MUS05_HH1-2024-03-29.mseed"  # HH1 mapped as vertical (trZ)
#file_p  = base_path / "XJ_MUS05_HDH-2024-03-29.mseed"  # Pressure trace
###############################################################
file_h1 = base_path / "XJ_MUS03_HH2-2023-09-22.mseed"
file_h2 = base_path / "XJ_MUS03_HH3-2023-09-22.mseed"
file_z  = base_path / "XJ_MUS03_HH1-2023-09-22.mseed"  # HH1 mapped as vertical (trZ)
file_p  = base_path / "XJ_MUS03_HDH-2023-09-22.mseed"  # Pressure trace
##########################################################################

################################################################
#file_h1 = base_path / "XJ_MUS05_HH2-2024-01-13.mseed"
#file_h2 = base_path / "XJ_MUS05_HH3-2024-01-13.mseed"
#file_z  = base_path / "XJ_MUS05_HH1-2024-01-13.mseed"
#file_p  = base_path / "XJ_MUS05_HDH-2024-01-13.mseed"


# ============================================================
# 4. Read miniSEED files
# ============================================================

tr1 = read(str(file_h1))[0]
tr2 = read(str(file_h2))[0]
trZ = read(str(file_z))[0]
trP = read(str(file_p))[0]


# ============================================================
# 5. Initialize DayNoise
# ============================================================

daynoise = DayNoise(

    tr1=tr1,
    tr2=tr2,
    trZ=trZ,
    trP=trP,

    # 2-hour analysis windows
    window=7200.0,

    # 30% overlap
    overlap=0.3,

    # Station identifier
    #key="XJ.MUS05"
    key="XJ.MUS03"
)


# ============================================================
# 6. Quality control
# ============================================================

daynoise.QC_daily_spectra(

    # Frequency passband in Hz
    #
    # 0.004 Hz = 250 s
    # 0.200 Hz =   5 s
    pd=[0.004, 0.2],

    # QC tolerance
    tol=1.5,

    fig_QC=False
)


# ============================================================
# 7. Calculate average spectra and tilt
# ============================================================

daynoise.average_daily_spectra(

    calc_rotation=True,

    # Tilt-estimation frequencies MUST remain in Hz
    #
    # 0.005 Hz = 200 s
    # 0.035 Hz = 28.57 s
    #tiltfreqs=[0.005, 0.035],
    tiltfreqs=[0.005, 0.025],

    # Generate tilt diagnostic figure
    fig_tilt=True,

    # No separate average-spectrum figure
    fig_average=False,

    # Do not save yet
    save=None,

    form="pdf"
)


# ============================================================
# 8. Get figure produced by DayNoise
# ============================================================

#fig = plt.gcf()

##############wamba ############################
fig = plt.gcf()
# Change H1 label to H2
for text in fig.findobj(match=plt.Text):

    current_text = text.get_text()

    if "CW from H1" in current_text:

        text.set_text(
            current_text.replace(
                "CW from H1",
                "CW from H2"
            )
        )
#######################################################
# ============================================================
# 9. Define period range and ticks
# ============================================================
#
# Original QC frequency range:
#
#     0.004 - 0.2 Hz
#
# corresponds to:
#
#     250 - 5 s
#
# We display:
#
#     10, 50, 100, 200 s
#
# ============================================================

period_min = 5.0
period_max = 250.0

period_ticks = np.array([
    10.0,
    50.0,
    #100.0,
    200.0
])


# ============================================================
# 10. Convert spectral plots from frequency to period
# ============================================================

for ax in fig.axes:

    # --------------------------------------------------------
    # Identify frequency/spectral axes
    # --------------------------------------------------------

    xlabel = ax.get_xlabel().lower()

    xmin, xmax = ax.get_xlim()

    xlower = min(xmin, xmax)
    xupper = max(xmin, xmax)

    is_spectral_axis = False


    # Explicit frequency-axis labels
    if (
        "frequency" in xlabel
        or "freq" in xlabel
        or "hz" in xlabel
    ):

        is_spectral_axis = True


    # Some DayNoise panels might not explicitly say
    # "Frequency", so identify relevant log-frequency axes.
    elif ax.get_xscale() == "log":

        if (
            xlower < 0.2
            and
            xupper > 0.004
        ):

            is_spectral_axis = True


    # Skip non-spectral axes
    if not is_spectral_axis:
        continue


    # ========================================================
    # 11. Convert each spectral line: frequency -> period
    # ========================================================

    for line in ax.get_lines():

        try:

            frequency = np.asarray(
                line.get_xdata(),
                dtype=float
            )

            ydata = np.asarray(
                line.get_ydata()
            )


            # Make sure x and y have the same length
            if len(frequency) != len(ydata):
                continue


            # ------------------------------------------------
            # Remove invalid frequencies
            #
            # T = 1/f cannot be calculated for f <= 0.
            # ------------------------------------------------

            valid = (
                np.isfinite(frequency)
                &
                (frequency > 0)
            )


            if not np.any(valid):
                continue


            frequency_valid = frequency[valid]

            y_valid = ydata[valid]


            # ------------------------------------------------
            # Convert frequency to period
            # ------------------------------------------------

            period = 1.0 / frequency_valid


            # ------------------------------------------------
            # Sort data by increasing period
            #
            # This results in:
            #
            # short period ------> long period
            # ------------------------------------------------

            order = np.argsort(period)

            period = period[order]

            y_valid = y_valid[order]


            # ------------------------------------------------
            # Replace frequency x-data with period x-data
            # ------------------------------------------------

            line.set_xdata(period)

            line.set_ydata(y_valid)


        except (TypeError, ValueError, IndexError):

            continue


    # ========================================================
    # 12. Set TRUE logarithmic PERIOD axis
    # ========================================================

    ax.set_xscale("log")

    #########change labels ####################
    ## ========================================================
    # 13. Set period limits
    # ========================================================
    #
    # Short period on left
    # Long period on right
    #
    #        5 s ------------> 250 s
    #
    # ========================================================

    ax.set_xlim(
        period_min,
        period_max
    )


    # ========================================================
    # 14. Define period ticks
    # ========================================================

    ax.set_xticks(
        period_ticks
    )


    # ========================================================
    # 15. Generate period tick labels automatically
    # ========================================================

    ax.set_xticklabels(

        [
            f"{period:g}"
            for period in period_ticks
        ],

        fontsize=9
    )


    # ========================================================
    # 16. Remove minor tick labels
    # ========================================================

    ax.tick_params(

        axis="x",

        which="minor",

        labelbottom=False
    )


    # ========================================================
    # 17. Change x-axis title
    # ========================================================

    ax.set_xlabel(

        "Period (s)",

        fontsize=10,

        labelpad=5
    )


# ============================================================
# 18. Adjust figure layout
# ============================================================

fig.tight_layout()


# ============================================================
# 19. Define output file
# ============================================================

output_file = (

    output_dir

    /

    #"XJ_MUS05_2024-01-13_tilt_angle_period_log.pdf"
    "XJ_MUS05_2023-09-22_tilt_angle_period_log.pdf"
    #"XJ_MUS05_2024-03-29_tilt_angle_period_log.pdf"
)


# ============================================================
# 20. Save figure
# ============================================================

fig.savefig(

    output_file,

    dpi=300,

    bbox_inches="tight"
)


# ============================================================
# 21. Close figure
# ============================================================

plt.close(fig)


# ============================================================
# 22. Print processing information
# ============================================================

print()
print("=============================================")
print("DayNoise / ATaCR processing completed")
print("=============================================")

print(
    f"Station: {daynoise.key}"
)

print(
    "X-axis: Period (s), logarithmic scale"
)

print(
    f"Period limits: "
    f"{period_min:g} - {period_max:g} s"
)

print(
    "Period ticks: "
    + ", ".join(
        f"{p:g}"
        for p in period_ticks
    )
    + " s"
)

print(
    f"Figure saved to: "
    f"{output_file.resolve()}"
)


# ============================================================
# 23. Print tilt direction if available
# ============================================================

if hasattr(daynoise, "rotation"):

    if hasattr(daynoise.rotation, "angle"):

        print(
            f"Tilt direction: "
            f"{daynoise.rotation.angle:.2f} degrees"
        )


print("=============================================")
print()
