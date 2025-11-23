import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import warnings
import numpy as np
import os
import sys
import time
from mpi4py import MPI
import pandas as pd
from Patrick.src.Utilities import set_inital_values,get_repo_path


class SimManager:
    def __init__(self, parameters, mesh_path, save_path=None,
                 parallel=False, runname: str = "initial_run", plot_only_run: bool = True, replace: bool = False):
        """
         Manages the configuration, setup, and execution of biochemical simulations.

         The `SimManager` class provides a high-level interface for initializing,
         configuring, and running simulations based on user-defined parameters.
         It is responsible for loading simulation models, setting up the required environment,
         and executing the simulation runs while saving the results in an HDF5 file.

         Args:
             parameters (dict): Dictionary containing the simulation parameters like
                                diffusion constants, reaction rates, etc.
             species_names (list): List of molecular species names involved in the simulation.
             mesh_path (str): File path to the geometry mesh used in the simulation.
             save_path (str): Path where simulation results will be stored.
             parallel (bool): Whether to run the simulation in parallel mode.
             runname (str): Identifier for the simulation run.
             plot_only_run (bool): If True, sets up the simulation for interactive
                     plotting only, without saving data. Default is False.
             replace (bool): If True, replaces existing simulation results with a new one.

         Attributes:
             parameters (dict): Dictionary containing the simulation parameters.
             species_names (list): List of molecular species names.
             endtime (int): Total simulation runtime, in seconds.
             mesh_path (str): File path to the geometry mesh.
             save_path (str): Path where simulation results will be stored.
             simulation (object): The initialized simulation instance.
             result_selector: Object for selecting and accessing simulation results.
             mesh: The loaded mesh geometry.
             cell_tets: Tetrahedral elements representing the cell.
             nuc_tets: Tetrahedral elements representing the nucleus.
             parallel (bool): Whether parallel mode is enabled.
             runname (str): Identifier for the simulation run.

        Methods:
            __init__(parameters, species_names, endt, mesh_path, save_file):
                Initializes the simulation manager and sets up the environment.
            _setup_environment():
                Ensures the necessary directory structure exists for saving simulation results.
            load_model(type):
                Loads a simulation model based on the provided type ("small" or "large").
            run(run_id=0):
                Executes the simulation and saves the results to an HDF5 file.
        """

        # Dynamically resolve the base repository path
        self.base_path = get_repo_path()

        # Resolve paths dynamically
        self.parameters = parameters
        self.species_names = None
        self.model_dataframe = None
        self.endtime = self.parameters["endtime"]
        self.mesh_path = mesh_path# Resolve mesh path
        #self.mesh_path = os.path.join(self.base_path, mesh_path)  # Resolve mesh path
        #self.save_path = os.path.join(self.base_path, save_path) if save_path else os.path.join(self.base_path, "saved_objects/initial_run/initial_run.h5")
        self.save_path = save_path if save_path else f"{self.base_path}saved_objects/initial_run/initial_run.h5"
        self.plot_only_run = plot_only_run
        self.replace = replace
        self.simulation = None
        self.result_selector = None
        self.mesh = None
        self.cell_tets = None
        self.nuc_tets = None
        self.parallel = parallel
        self.runname = runname

        self._setup_environment(self.save_path)
        self._load_model_dataframe()

    def _load_model_dataframe(self):
        """
        Load and clean the model dataframe dynamically based on the parameters.
        """
        from src.Utilities import dataframe_cleanup
        df_path = f"{self.base_path}{self.parameters["big_model_mini_sph_df_path"]}"
        df = pd.read_excel(df_path)
        df = dataframe_cleanup(df, ["Species"])
        self.model_dataframe = df
        self.species_names = df.Species.values

    def _setup_environment(self, path):
        """
        Set up directories and environment.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        save_dir = os.path.dirname(path)  # Get the parent directory of the path

        if rank == 0: # makes sure that only one process does this
            if not os.path.exists(save_dir):  # Check if the directory exists
                os.makedirs(save_dir, exist_ok=True)  # Create it if it doesn't exist

            if os.path.isfile(path + ".h5") and self.replace and not self.plot_only_run:
                try:
                    os.remove((path + ".h5"))
                except FileNotFoundError:
                    pass  # If the file was already removed by another process, ignore it


    def load_model(self, type, mesh_scale=1):
        """
        Load a simulation model based on the specified type.

        Args:
            type (str): Type of model to load, currently supports "small" with
                       "large" planned for future implementation.

        Raises:
            UserWarning: If the specified model type is not implemented.
        """
        if type == "small":
            from src.Model_small import create_model
            self.simulation, self.result_selector, self.mesh = create_model(self.parameters,
                                                                            self.species_names,
                                                                            self.mesh_path,
                                                                            mesh_scale,
                                                                            self.plot_only_run)
            # self.model_data = pd.read_excel("/home/pb/steps_cell_signaling/Patrick/data_small_model.xlsx")
        elif type == "large":
            from src.Model_expanded_mini_sph_new import create_model
            # warnings.warn("The 'large' model type is not yet implemented", UserWarning)
            self.simulation, self.result_selector, self.mesh = create_model(self.model_dataframe,
                                                                            self.parameters,
                                                                            self.species_names,
                                                                            self.mesh_path,
                                                                            self.plot_only_run)
        else:
            warnings.warn(f"The '{type}' model type is not yet implemented", UserWarning)

    def run(self, replicats):
        """
        Run the simulation with the initialized parameters.

        This method either saves results to an HDF5 file or provides interactive
        plotting during the simulation run, based on the plot_only_run setting.

        Args:
            replicats (int): How many replicats of this sim to run. Default is 1.

        Notes:
            - When plot_only_run is False, results are saved to the specified HDF5 file.
            - When plot_only_run is True, an interactive plotting session is launched.
            - Initial molecular counts are set for EGF, EGFR, and GAP species regardless of mode.
        """
        # this is a regex way of finding out if a file with a certain fileformat already exists. Right now its hardcoded
        # that we always use .h5. Would be nice to make it dynamic
        # p = "../file.h5"
        # fileformat = ["h5","csv"]
        #
        # import re
        # m = re.match(fr".+\.(?:{'|'.join(fileformat)})$")
        # if m is not None:

        if self.plot_only_run == False:
            # if replicat_id:
            #     replicat_save_path = self.save_path + f"/replicat_{replicat_id}/output"
            #     self._setup_environment(replicat_save_path)
            #     checked_save_path = replicat_save_path
            # else:
            checked_save_path = self.save_path

            sum_of_runtimes = 0
            for i in range(replicats):
                options = dict(compression="gzip", compression_opts=5) # compress the output files to save space. Larger opts = more compression
                with stsave.XDMFHandler(checked_save_path, hdf5DatasetKwArgs=options) as hdf:
                    self.simulation.toDB(hdf, uid = self.runname)
                    self.simulation.newRun()
                    
                    set_inital_values(self, factor = self.initial_factor) # TODO: MAYBE DONT HARDCODE BRO

                    # self.simulation.exo.EGF.Count = self.parameters["EGF_0"]
                    # self.simulation.cell_surface.EGFR.Count = self.parameters["EGFR_0"]
                    # self.simulation.cyt.GAP.Count = self.parameters["GAP_0"]
                    # self.simulation.cyt.ERK.Count = self.parameters["ERK_0"]
                    # self.simulation.cyt.P3.Count = self.parameters["P3_0"]

                    self.simulation.nuc_mem.ERKp.DiffusionActive = True

                    start_time = time.time()
                    self.simulation.run(self.endtime)
                    end_time = time.time()

                    print(f"Run completed in {end_time - start_time:.2f} seconds.")
                    sum_of_runtimes += (end_time - start_time)
            print(f"All runs completed in {sum_of_runtimes:.2f} seconds")
        else:
            comm = MPI.COMM_WORLD
            assert comm.Get_size() == 1, "A plot_only_run only works in serial due to limitations by STEPS. Rerun either without mpirun or with mpirun -n 1 ..."

            from src.InteractivePlotting import interactive_plots
            self.simulation.newRun()
            self.simulation.exo.EGF.Count = self.parameters["EGF_0"]
            self.simulation.cell_surface.EGFR.Count = self.parameters["EGFR_0"]
            self.simulation.cyt.GAP.Count = self.parameters["GAP_0"]
            self.simulation.cyt.ERK.Count = self.parameters["ERK_0"]
            self.simulation.cyt.P3.Count = self.parameters["P3_0"]

            self.simulation.nuc_mem.ERKp.DiffusionActive = True

            self.result_selector = stsave.ResultSelector(self.simulation)
            SimControl = interactive_plots(self)
            SimControl.run()

