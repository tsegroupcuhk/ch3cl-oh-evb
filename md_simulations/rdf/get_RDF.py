import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

TINY_SIZE = 8
SMALL_SIZE = 8
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=TINY_SIZE, titleweight='bold')  # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=TINY_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=TINY_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=TINY_SIZE-1.0)  # legend fontsize


# Load your trajectory and topology file
# Replace 'topology_file' and 'trajectory_file' with your actual file paths
u = mda.Universe("sys.psf","rdf.dcd")

# Select atom groups (example: oxygen and hydrogen atoms in water)
oxygen = u.select_atoms("type OOHR")
hydrogen = u.select_atoms("type ODW")

# Create RDF object
rdf = InterRDF(oxygen, hydrogen, nbins=100, range=(0.0, 10.0))

# Run RDF analysis using every frame
rdf.run(start=0, stop=50_000,step=1) # total 50,000 frames (50 ns)

# Find peaks in the RDF
peaks, _ = find_peaks(rdf.rdf)

# Check if any peaks are found
if len(peaks) > 0:
    first_peak_index = peaks[0]
    first_peak_position = rdf.bins[first_peak_index]
    first_peak_height = rdf.rdf[first_peak_index]
    print(f"First RDF peak at r = {first_peak_position:.2f} Å, height = {first_peak_height:.2f}")
else:
    print("No peaks found in RDF.")

def moving_average(data, window_size):
    # Create a window (kernel) of ones, normalized by the window size
    window = np.ones(window_size) / window_size
    # Convolve the data with the window
    smoothed_data = np.convolve(data, window, mode='same')
    return smoothed_data
window_size = 1
y = moving_average(rdf.rdf, window_size)

# Optional: Plot RDF and mark the first peak
fig = plt.figure(figsize=(3.54,2.54), dpi=600)
plt.plot(rdf.bins, y, label="O-O* RDF",color='black',linewidth=1)
if len(peaks) > 0:
    plt.plot(first_peak_position, first_peak_height, "ro",markersize=3)
    plt.annotate(f'({first_peak_position:.2f}, {first_peak_height:.2f})', (first_peak_position, first_peak_height), textcoords="offset points", xytext=(8, 0))

plt.xlabel(r"$r (\mathrm{Å})$")
plt.ylabel(r"$g(r)$")
#plt.legend(frameon=False)
plt.xlim(0, 10)
plt.ylim(-0.2, 2.5)
plt.savefig("O-Ostar-RDF.png", bbox_inches='tight')

# Save RDF data to a text file for Origin
with open("RDF_data.txt", "w") as f:
    f.write("# r (Å)\tg(r)\n")  # Header
    for bin_position, rdf_value in zip(rdf.bins, rdf.rdf):
        f.write(f"{bin_position:.6f}\t{rdf_value:.6f}\n")
