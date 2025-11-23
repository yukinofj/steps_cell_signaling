import sys
from src.SimManager import SimManager
import logging
from parameters import p
from Patrick.src.Utilities import get_repo_path
import numpy as np
import steps.saving as stsave
import os
# cd ~/code/steps_cell_signaling/; PYTHONPATH=$PWD python3 Patrick/run3.py


ellipsoidity = 1.0      # ellipsoidity: 0=sphere, 1=ellipsoid
n_reps = 100             # number of repetions
runnr = 9
all_results = []
run = "testrun"         # run name
parameters = list(p.keys())
dt = p[parameters[0]]
t_end = p[parameters[1]]
self.initial_factor = 1.0

file_name = f"{run}{runnr}_E{ellipsoidity}_N{n_reps}_dt{dt}_tend{t_end}"

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logging.info("Starting the script...")

try:
    # Retrieve the base repository path
    logging.info("Retrieving the repository path...")
    base_path = get_repo_path()
    logging.info(f"Repository path retrieved: {base_path}")

    save_path = f"{base_path}Patrick/saved_objects/run{runnr}/{file_name}"

    # Initialize the simulation manager
    logging.info("Initializing the simulation manager...")
    sm = SimManager(parameters=p,
        mesh_path = f"{base_path}Patrick/meshes_ellipsoidity/ellipsoidity_{ellipsoidity}.inp",
        save_path= save_path,  # Full file path without the .h5 suffix
        parallel = True, # no effect so far
        runname = f"{run}",
        plot_only_run = False,
        replace = True) # whether an already existing file should be overwritten or not. Might throw an error if there is an already existing one and this is set to false.
    logging.info("Simulation manager initialized.")
    # Load the model
    logging.info("Loading the model...")
    sm.load_model(type="small", mesh_scale=1)
    logging.info("Model loaded successfully.")

    # Run the simulation
    logging.info("Running the simulation...")
    sm.run(replicats = n_reps)
    logging.info("Simulation completed successfully.")
except Exception as e:
    logging.error(f"An error occurred: {e}", exc_info=True)