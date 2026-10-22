from pathlib import Path
import matplotlib.pyplot as plt
from obspy import read
from obstools.atacr import DayNoise
import numpy as np
# 1. Define paths and output setup
#base_path = Path("LucerneData/muota-SingSta-05/mseed/2024/03/29")
#base_path = Path("LucerneData/muota-SingSta-03/mseed/2023/09/22")
base_path = Path("LucerneData/muota-SingSta-05/mseed/2024/01/13")
output_dir = Path("./figures")
output_dir.mkdir(parents=True, exist_ok=True)

## File names
#file_h1 = base_path / "XJ_MUS05_HH1-2024-03-29.mseed"
#file_h2 = base_path / "XJ_MUS05_HH2-2024-03-29.mseed"
#file_z  = base_path / "XJ_MUS05_HH3-2024-03-29.mseed"  # HH3 mapped as vertical (trZ)
#file_p  = base_path / "XJ_MUS05_HDH-2024-03-29.mseed"  # Pressure trace
#################################################################################

#file_h1 = base_path / "XJ_MUS03_HH1-2023-09-22.mseed"
#file_h2 = base_path / "XJ_MUS03_HH2-2023-09-22.mseed"
#file_z  = base_path / "XJ_MUS03_HH3-2023-09-22.mseed"  # HH3 mapped as vertical (trZ)
#file_p  = base_path / "XJ_MUS03_HDH-2023-09-22.mseed"  # Pressure trace
######################################################################################
file_h1 = base_path / "XJ_MUS05_HH1-2024-01-13.mseed"
file_h2 = base_path / "XJ_MUS05_HH2-2024-01-13.mseed"
file_z  = base_path / "XJ_MUS05_HH3-2024-01-13.mseed"  # HH3 mapped as vertical (trZ)
file_p  = base_path / "XJ_MUS05_HDH-2024-01-13.mseed"  # Pressure trace

# 2. Read miniSEED files into ObsPy Trace objects
tr1 = read(str(file_h1))[0]
tr2 = read(str(file_h2))[0]
trZ = read(str(file_z))[0]
trP = read(str(file_p))[0]

# 3. Initialize DayNoise object
daynoise = DayNoise(
    tr1=tr1,
    tr2=tr2,
    trZ=trZ,
    trP=trP,
    window=7200.0,  # 2-hour windows (7200 s)
    overlap=0.3,    # 30% overlap
    key="XJ.MUS05"
    #key="XJ.MUS03"
)

# 4. Perform Quality Control on daily noise spectra
daynoise.QC_daily_spectra(
    pd=[0.004, 0.2], # Passband corners (Hz)
    tol=1.5,         # Tolerance threshold
    fig_QC=False
)

# 5. Calculate tilt angle (disable default plot generation)
daynoise.average_daily_spectra(
    calc_rotation=True,
    tiltfreqs=[0.005, 0.035], # Passband (Hz) for tilt estimation
    #tiltfreqs=[0.005, 0.02], # Passband (Hz) for tilt estimation
    fig_tilt=False,
    fig_average=False
)

# 6. Extract rotation data and plot Coherence vs. Tilt Angle

# 6. Extract rotation data and plot Coherence vs. Tilt Angle
if hasattr(daynoise, 'rotation') and daynoise.rotation is not None:
    rot = daynoise.rotation

    # Safe extraction of coherence array
    if hasattr(rot, 'coh'):
        coherence = rot.coh
    elif isinstance(rot, dict) and 'coh' in rot:
        coherence = rot['coh']
    else:
        raise AttributeError("Could not find coherence array in daynoise.rotation")

    # Safe extraction of optimal tilt angle
    if hasattr(rot, 'tilt'):
        optimal_angle = rot.tilt
    elif hasattr(rot, 'angle'):
        optimal_angle = rot.angle
    elif isinstance(rot, dict):
        optimal_angle = rot.get('tilt', rot.get('angle', None))
    else:
        optimal_angle = None

    # Reconstruct the search angles (0 to 360 degrees) matching the coherence array length
    angles = np.linspace(0, 360, len(coherence), endpoint=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(angles, coherence, color='k', lw=2, label='Coherence')

    # Highlight optimal tilt angle if found
    if optimal_angle is not None:
        ax.axvline(optimal_angle, color='r', linestyle='--', label=f'Optimal Tilt: {optimal_angle:.2f}°')

    ax.set_xlabel('Tilt Angle (Degrees from H1)')
    ax.set_ylabel('Coherence')
    ax.set_title(f'Coherence vs. Tilt Angle — Station {daynoise.key}')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='best')

    # Save figure
    #output_file = output_dir / "XJ_MUS05_2024-03-29_coherence_vs_tilt.pdf"
    output_file = output_dir / "XJ_MUS05_2024-01-13_coherence_vs_tilt.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Station: {daynoise.key}")
    if optimal_angle is not None:
        print(f"Tilt Direction (Azimuth from H1): {optimal_angle:.2f}°")
    print(f"Saved custom coherence plot to: {output_file.resolve()}")
