import steps.interface
import steps.saving as stsave
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
#import matplotlib.cm as cm
import seaborn as sns
import h5py
import numpy as np
import re
import math
from matplotlib import rc
from parameters import p

rc("text", usetex=True)
#rc("font", family="sans-serif")
#rc("text.latex", preamble=r"\usepackage{sfmath\renewcommand{\rmdefault}{cmss}}")

# parameters
name = "EGF"
key = "EGF_initial_conc"
runnr =9
n_reps = 50
ellipsoidity = 0.0
base_c = p.get(f"{key}", 1.0)
c = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
dt = p["time step"]
t_end = p["endtime"]

input_dir = f"/home/yukinofj/code/steps_cell_signaling/Patrick/saved_objects/parameter_scan/run{runnr}"
cmap = plt.colormaps["cool"]
norm = mcolors.Normalize(vmin = min(c), vmax = max(c))

hdf_path_1 = f"{input_dir}/PSrun{runnr}_{name}{c[0]}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}" # PSrun{runnr}_{name}{fraction}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}
hdf = stsave.HDF5Handler(hdf_path_1)
results = hdf[name].results
full_labels = [x.labels for x in results]
species_names = [re.search(r'\.(.*?)\.', s[0]).group(1) for s in full_labels if re.search(r'\.(.*?)\.', s[0])]
species_latex = {
    "EGF":              r"EGF",
    "EGF_EGFR":         r"EGF-EGFR",
    "EGF_EGFR2":        r"(EGF-EGFR)$^2$",
    "EGF_EGFRp2":       r"(EGF-EGFR$^p$)$^2$",
    "EGF_EGFRp2_GAP":   r"EGF-EGFR$^p$)$^2$-GAP",
    "ERK":              r"ERK",
    "ERKp":             r"ERK$^p$",
    "ERKp":             r"ERK$^p$",
    "ERKpp":            r"ERK$^{pp}$"
}

# Plot all results
n_species = len(species_names)
grid_size = math.ceil(math.sqrt(n_species))  # Prefer a square layout
n_rows, n_cols = grid_size, math.ceil(n_species / grid_size)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6, 4), constrained_layout=True)
axes = axes.flatten()  # Flatten in case of 2D array
for ax in axes.flat:
    ax.tick_params(axis='both', which='both', direction='in')

for idx, (res, species_name) in enumerate(zip(results, species_names)):
    ax = axes[idx]
    for c_v in c:
        hdf_path = f"{input_dir}/PSrun{runnr}_{name}{c_v}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"
        hdf = stsave.HDF5Handler(hdf_path)
        results = hdf[name].results
        res = results[idx]
        data_start_mean = np.mean(res.data[:,0,0])
        #print(f"c_v={c_v}, initial mean data )= {data_start_mean}")
        mean_data = np.mean(res.data[:, :, 0], axis=0)
        std_data = np.std(res.data[:, :, 0], axis=0)
        ax.plot(res.time[0], mean_data, color=cmap(norm(float(c_v))), label=f"{c_v}")
        ax.fill_between(res.time[0], mean_data - std_data, mean_data + std_data, color=cmap(norm(float(c_v))), alpha=0.25, edgecolor="none")
    if idx >= n_species - n_cols:
        ax.set_xlabel('Time [s]')
    label = species_latex.get(species_name, species_name)
    ax.set_ylabel(label)
# Hide any unused subplots

for i in range(idx + 1, len(axes)):
    fig.delaxes(axes[i])


sm = plt.cm.ScalarMappable (cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
cbar.set_label(f"{name}", rotation=270, labelpad=15)
output = f"/home/yukinofj/code/steps_cell_signaling/Patrick/figures/plot_PS/plot_PSrun{runnr}_{name}{c_v}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}.pdf"
plt.savefig(output, bbox_inches="tight", transparent=True)
print(f"Plot saved to {output}.")

