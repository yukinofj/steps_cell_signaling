import os
from os.path import join
from src.SimManager import SimManager
from parameters import p
import numpy as np
import pandas as pd
import steps.API_2.geom as stgeom
# instead of running over ellipsoidity meshes, we can loop over concentrations, for example
# read base count from excel file; modify value for each scan iteration; update simulations initial conditions before running
# cd ~/code/steps_cell_signaling; PYTHONPATH=$PWD python3 Patrick/run_parameter_scan.py


# parameters
name = "EGF"
key = "EGF_initial_conc"
runnr = 12
n_reps = 1
ellipsoidity = 1.0
base_c = p.get(f"{key}", 1.0)
fractions = [0.1]#[0.8, 0.9, 0.95, 0.975, 0.99, 1.0]
dt = p["time step"]
t_end = p["endtime"]
home_dir = "/home/yukinofj/code/steps_cell_signaling/"

# read base count
xls_path = os.path.join(home_dir, p["big_model_mini_sph_df_path"])
df = pd.read_excel(xls_path)
c_row = df[df["Species"].str.strip() == f"{name}'"]
base_count = c_row["exo init count"].values[0]
print(f"Base count from excel: {base_count}")

# save dir
save_dir = os.path.join(home_dir, f"Patrick/saved_objects/parameter_scan_{name}/run{runnr}")
os.makedirs(save_dir, exist_ok=True)
print(f"Directory {save_dir} created.")
save_path = f"{save_dir}/PSrun{runnr}_{name}{fractions}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"

def set_initial_values(sim_manager, factor):
    # read excel
    df_path = os.path.join(sim_manager.base_path, sim_manager.parameters["big_model_mini_sph_df_path"])
    df = pd.read_excel(df_path)
    df["Species"] = df["Species"].str.strip() # strip whitespace

    # construct compartment-to-column mapping (names ending in "init count")
    # build dictionary with name key
    init_cols = df.filter(like="init count").columns
    comp_map = {col.split("init count")[0].strip(): col for col in init_cols}
    compartments = {}
    for key, value in sim_manager.simulation._children.items():
        if isinstance(value, stgeom._TetCompartment) or isinstance(value, stgeom._TetPatch):
            compartments[key] = value

    # set initial values for each species in each compartment
    for comp_name, col_name in comp_map.items():
        for species in sim_manager.species_names:
            try:
                initial_value = df.loc[df["Species"] == species, col_name].values[0]
            except IndexError:
                continue
            initial_value = 0 if np.isnan(initial_value) else initial_value
            value = initial_value * factor
            try:
                comp_obj = compartments[comp_name]
                getattr(comp_obj, species).Count = value
                print(comp_name, species, value)
            except steps.API_2.sim.SolverCallError:
                pass 

# run sim for each parameter in "fractions"
for frac in fractions:
    file_name = f"PSrun{runnr}_{name}{frac}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"
    save_path = os.path.join(save_dir, file_name)

    sm = SimManager(parameters=p,
                    mesh_path=f"/home/yukinofj/code/steps_cell_signaling/Patrick/meshes_ellipsoidity/ellipsoidity_{ellipsoidity}.inp",
                    save_path=save_path,
                    parallel=True,
                    runname=name,
                    plot_only_run=False,
                    replace=True)
    sm.load_model(type="small", mesh_scale=1)
    sm.initial_factor = frac
    sm.run(replicats = n_reps)
    print(f"Scan for {frac} completed and saved in {save_path}.")