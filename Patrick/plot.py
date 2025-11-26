import steps.interface
import steps.API_2.saving as stsave
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
#import matplotlib.cm as cm
import seaborn as sns
import h5py
import numpy as np
import re
import math
from matplotlib import rc
from parameters import p
import argparse

# terminal input
# cd /home/yukinofj/code/steps_cell_signaling/Patrick; python3 plot.py 9

def format_param(val, decimals=4):
    return float(f"{val:.{decimals}g}")


pars = argparse.ArgumentParser()
pars.add_argument("runnr", type=int)
args = pars.parse_args()
print("runnr:", args.runnr)

runnr = args.runnr
run = p["run"]
n_reps = p["n_reps"]
ellip_value = [0.0, 1.0] 
dt = p["time step"]
t_end = p["endtime"]

home_dir = f"/home/yukinofj/code/steps_cell_signaling/Patrick/saved_objects/run{runnr}"
fig_dir = "/home/yukinofj/code/steps_cell_signaling/Patrick/figures/plot"

# latex formatting
rc("text", usetex=True)
species_latex = {
    "EGF":              r"EGF",
    "EGF_EGFR":         r"EGF-EGFR",
    "EGF_EGFR2":        r"(EGF-EGFR)$^2$",
    "EGF_EGFRp2":       r"(EGF-EGFR$^p$)$^2$",
    "EGF_EGFRp2_GAP":   r"EGF-EGFR$^p$)$^2$-GAP",
    "ERK":              r"ERK",
    "ERKp":             r"ERK$^p$",
    "ERKpp":            r"ERK$^{pp}$"
}

# extract the species names for the result_selector label via regex
def get_species_names(group):
    #file_name = f"{run}{runnr}_ellip{E}_reps{n_reps}_dt{dt}_tend{t_end}"
    #hdf_path = f"{home_dir}/{file_name}"
    #hdf = stsave.HDF5Handler(hdf_path)
    #group = hdf["testrun"]
    results = group.results
    labels = [r.labels[0] for r in results]
    species = [re.search(r'\.(.*?)\.', lbl).group(1) for lbl in labels if re.search(r'\.(.*?)\.', lbl)]
    return species

def load_ellip_file(E):
    file_name = f"{run}{runnr}_E{E}_N{n_reps}_dt{dt}_tend{t_end}"
    hdf_path = f"{home_dir}/{file_name}"
    hdf = stsave.HDF5Handler(hdf_path)
    group = hdf["testrun"]
    species_names = get_species_names(group) # hdf
    results = group.results
    means = []
    stds = []
    times = None

    for res in results:
        data = np.array(res.data)[:,:,0]
        if times is None:
            times = res.time[0]
        means.append(np.mean(data, axis = 0))
        stds.append(np.std(data, axis = 0))
    return species_names, times, np.array(means), np.array(stds)

species_names, _, _, _ = load_ellip_file(ellip_value[0]) # for species list
n_species = len(species_names)
print("Species names: ", species_names)

# subplot grid
grid_size = math.ceil(math.sqrt(n_species))
n_rows, n_cols = grid_size, math.ceil(n_species / grid_size)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6, 4), constrained_layout=True)
axes = axes.flatten()
for ax in axes.flat:
    ax.tick_params(axis='both', which='both', direction='in')
colors = ["deepskyblue", "deeppink", "blueviolet"]
# loop through species
for idx, species in enumerate(species_names):
    ax = axes[idx]
    # loop over ellip
    for i, ellip in enumerate(ellip_value):
        _, times, mean, std = load_ellip_file(ellip)
        mean_data = mean[idx]
        std_data = std[idx]
        color = colors[i]
        ax.plot(times, mean_data, label = f"$E ={ellip}$", color=color)
        lower_bound = np.clip(mean_data - std_data, 0, None)
        ax.fill_between(times, lower_bound, mean_data + std_data, color=color, alpha=0.25, edgecolor="none")
    label = species_latex.get(species, species)
    ax.set_ylabel(label)
    if idx >= n_species - n_cols:
        ax.set_xlabel('Time [s]')
    ax.legend(fontsize= "small", frameon = False)
# Hide any unused subplots
for i in range(idx + 1, len(axes)):
    fig.delaxes(axes[i])

ellip_str = "_".join(str(x) for x in ellip_value)
output = f"{fig_dir}/plot_{run}{runnr}_E{ellip_str}_N{n_reps}_dt{dt}_tend{t_end}.pdf"
plt.savefig(output, bbox_inches="tight", transparent=True)
print(f"Plot saved to {output}.")