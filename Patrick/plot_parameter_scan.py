import steps.interface
import steps.saving as stsave
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sns
import h5py
import numpy as np
import re
import math
from matplotlib import rc
from parameters import p
import argparse

# terminal input
# cd /home/yukinofj/code/steps_cell_signaling/Patrick; python3 plot_parameter_scan.py 9 "linspace(0,1,6)"

def format_param(val, decimals=4):
    return float(f"{val:.{decimals}g}")


def parse_fractions(arg):
    arg = arg.strip()
    # linspace(...) or logspace(...)
    if arg.startswith("linspace") or arg.startswith("logspace"):
        inside = arg.split("(", 1)[1][:-1]
        start, stop, n_points = inside.split(",")
        start, stop, n_points = float(start), float(stop), int(float(n_points))
        if arg.startswith("linspace"):
            vals = np.linspace(start, stop, n_points)
        else:
            vals = np.logspace(start, stop, n_points)
        return [float(v) for v in vals]

    # comma-separated list: "0.8,0.9,0.95,0.975,0.99,1.0"
    return [float(v) for v in arg.split(",")]


pars = argparse.ArgumentParser()
pars.add_argument("runnr", type=int)
pars.add_argument("fractions", type=parse_fractions)
args = pars.parse_args()
print("runnr:", args.runnr)
print("fractions:", args.fractions)

runnr, fractions = (args.runnr, args.fractions)

name = p["parameter scan"] # EGF
key = "EGF_initial_conc"
n_reps = p["n_reps"]
ellipsoidity = p["ellipsoidity"]
dt = p["time step"]
t_end = p["endtime"]

base_EGF = p.get(key, 1.0)
init_vals = [base_EGF * f for f in fractions]

rc("text", usetex=True)
#rc("font", family="sans-serif")
#rc("text.latex", preamble=r"\usepackage{sfmath\renewcommand{\rmdefault}{cmss}}")
'''
# parameters
name = "EGF"
key = "EGF_initial_conc"
runnr = 9
n_reps = 50
ellipsoidity = 0.0
base_value = p.get(f"{key}", 1.0)
fractions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
'''

input_dir = f"/home/yukinofj/code/steps_cell_signaling/Patrick/saved_objects/parameter_scan_EGF/run{runnr}"
cmap = plt.colormaps["cool"]
norm = mcolors.Normalize(vmin = min(fractions), vmax = max(fractions))

hdf_path_1 = f"{input_dir}/PSrun{runnr}_{name}{format_param(fractions[0])}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}" # PSrun{runnr}_{name}{fraction}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}
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
grid_size = math.ceil(math.sqrt(n_species))
n_rows, n_cols = grid_size, math.ceil(n_species / grid_size)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6, 4), constrained_layout=True)
fig2, axes2 = plt.subplots(1, 2, figsize=(5, 2), constrained_layout=True)
axes, axes2 = axes.flatten(), axes2.flatten()
for ax in axes.flat:
    ax.tick_params(axis='both', which='both', direction='in')
for ax in axes2.flat:
    ax.tick_params(axis='both', which='both', direction='in')

for idx, (res, species_name) in enumerate(zip(results, species_names)):
    ax = axes[idx]
    for i, c_v in enumerate(fractions):
        hdf_path = f"{input_dir}/PSrun{runnr}_{name}{format_param(c_v)}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"
        hdf = stsave.HDF5Handler(hdf_path)
        results = hdf[name].results
        res = results[idx]
        data_start_mean = np.mean(res.data[:,0,0])
        mean_data = np.mean(res.data[:, :, 0], axis=0)
        std_data = np.std(res.data[:, :, 0], axis=0)
        ax.plot(res.time[0], mean_data, color=cmap(norm(float(c_v))), label=f"{c_v}")
        lower_bound = np.clip(mean_data - std_data, 0, None)
        ax.fill_between(res.time[0], lower_bound, mean_data + std_data, color=cmap(norm(float(c_v))), alpha=0.25, edgecolor="none")
    if idx >= n_species - n_cols:
        ax.set_xlabel('Time [s]')
    label = species_latex.get(species_name, species_name)
    ax.set_ylabel(label)
for i in range(idx + 1, len(axes)):
    fig.delaxes(axes[i])

erkpp_peaks = []
for i, c_v in enumerate(fractions):
    hdf_path = f"{input_dir}/PSrun{runnr}_{name}{format_param(c_v)}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"
    hdf = stsave.HDF5Handler(hdf_path)
    results = hdf[name].results
    color = cmap(norm(c_v))

    egf_raw = results[species_names.index(name)]
    egf_val = np.mean(egf_raw.data[:,:,0],axis=0)
    erkpp_raw = results[species_names.index("ERKpp")] # how "name" influences signal/ERKpp
    erkpp_val = np.mean(erkpp_raw.data[:,:,0],axis=0)
    peak_val = np.max(erkpp_val)
    erkpp_peaks.append(peak_val)

    axes2[1].semilogx(erkpp_val,egf_val, color=cmap(norm(float(c_v))), label=f"{c_v}")
    axes2[1].set_ylabel("EGF")
    axes2[1].set_xlabel(r"ERK$^{pp}$")

for i, peak in enumerate(erkpp_peaks):
    color = cmap(norm(fractions[i]))
    axes2[0].loglog(init_vals[i], peak, "o", color=color)
axes2[0].set_xlabel("EGF init con")
axes2[0].set_ylabel(r"max ERK$^{pp}$")



sm = cm.ScalarMappable (cmap=cmap, norm=norm)
sm2 = cm.ScalarMappable (cmap=cmap, norm=norm)
sm.set_array([])
sm2.set_array(init_vals)
cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
cbar2 = fig2.colorbar(sm2, ax=axes2, orientation="vertical", fraction=0.02, pad=0.04)
cbar.set_label(f"{name}", rotation=270, labelpad=15)
cbar2.set_label(f"{name} init conc", rotation=270, labelpad=15)
output = f"/home/yukinofj/code/steps_cell_signaling/Patrick/figures/plot_PS/signpath_PSrun{runnr}_{name}{format_param(c_v)}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}.pdf"
output2 = f"/home/yukinofj/code/steps_cell_signaling/Patrick/figures/plot_PS/erkpp_PSrun{runnr}_{name}{format_param(c_v)}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}.pdf"
fig.savefig(output, bbox_inches="tight", transparent=True)
fig2.savefig(output2, bbox_inches="tight", transparent=True)
print(f"Plot saved to {output}.")

