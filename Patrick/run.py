import sys
from src.SimManager import SimManager
import logging
from parameters import p
from Patrick.src.Utilities import get_repo_path

"""
Initialize and run the simulation for cell signaling pathways.

This script sets up a simulation manager instance with specified parameters and species,
then loads a small model and executes a single simulation run with parallel processing.

The simulation models EGF-EGFR signaling pathway dynamics within a spherical cell mesh.

If you run a "plot_only_run" the result_selectors need to be different whether you want to plot the results on the go or 
save them. Setting "plot_only_run = True" adjusts it accordingly. But: when you do a "plot_only_run" it means there is 
no data saved, so be aware.
"""
# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Starting the script...")

try:
    # Retrieve the base repository path
    logging.info("Retrieving the repository path...")
    base_path = get_repo_path()
    logging.info(f"Repository path retrieved: {base_path}")
   
    # # Get the replicat_id from the command-line argument
    # if len(sys.argv) < 2:
    #     replicat_id = None
    # else:
    #     replicat_id = int(sys.argv[1])
    # todo: implement relative
    # Initialize the simulation manager
    logging.info("Initializing the simulation manager...")
    sm = SimManager(parameters=p,
                    # mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/sph_1.0.inp",
                    # mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/elipsoid_4.5.inp",
                    mesh_path = f"{base_path}Patrick/meshes_ellipsoidity/ellipsoidity_0.0.inp",
                    #save_path ="/home/pb/steps_cell_signaling/Patrick/saved_objects/testing/test2",  #without the .h5 suffix, but full file path please
                    # mesh_path=f"{base_path}Patrick/meshes_ellipsoidity/ellipsoidity_0.8.inp",
                    save_path=f"{base_path}Patrick/saved_objects/ellipsoidity_runs/testrun1",  # Full file path without the .h5 suffix
                    parallel = True, # no effect so far
                    runname = "test",
                    plot_only_run = False,
                    replace = True) # whether an already existing file should be overwritten or not. Might throw an error if there is an already existing one and this is set to false.
    logging.info("Simulation manager initialized.")
    # Load the model
    logging.info("Loading the model...")
    sm.load_model(type="small", mesh_scale=1)
    logging.info("Model loaded successfully.")

    # Run the simulation
    logging.info("Running the simulation...")
    sm.run(replicats = 1)
    logging.info("Simulation completed successfully.")
except Exception as e:
    logging.error(f"An error occurred: {e}", exc_info=True)