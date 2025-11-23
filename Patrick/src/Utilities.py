from scipy.constants import N_A
import contextlib
import io
import os
import sys
import json
import pandas as pd
import numpy as np
import steps.geom as stgeom
import steps.interface
from os.path import abspath, dirname, join
import socket
import getpass

def molar_to_molecules(M, Volume):
    """
    Convert molar concentration to number of molecules.

    Args:
        M (float): Molar concentration in mol/L.
        Volume (float): Volume in cubic meters (m^3).

    Returns:
        float: Number of molecules in the given volume at the specified concentration.
    """
    return M * (N_A * Volume * 1e3)  # 1e3 bcs M = mol/l = 1000 mol/m^3


@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout

# cleans the dataframe, will probably be extended
def dataframe_cleanup(df, columns=["Species"]):
    for col in columns:
        if col == "Species":
            df[col] = df[col].str.replace("'", "", regex=False)
            df[col] = df[col].str.replace(" ", "", regex=False)
    return df



def set_inital_values(sim_manager, factor):
    '''
    Initializes species counts in compartments based on an Excel configuration file.

    This function reads initial species count data from an Excel file, cleans it,
    and maps the values to the appropriate compartments in the simulation.
    The initial values are scaled by a given factor before being assigned.

    Parameters
    ----------
    factor : float
        A multiplier applied to the initial species count values before assignment.

    Notes
    -----
    - Reads the data from an Excel file specified in `sim_manager.parameters["big_model_mini_sph_df_path"]`.
    - Cleans the dataframe and extracts species names.
    - Identifies compartments in the simulation where species should be initialized.
    - Assigns initial values from the dataframe to species in respective compartments.
    - If a value is NaN, it defaults to 0.
    - Any `SolverCallError` from STEPS API is caught and ignored.
    parameter scan implementieren, zb EGF (excel exo init count) 
    '''
    df_path = f"{sim_manager.base_path}{sim_manager.parameters["big_model_mini_sph_df_path"]}"
    df = pd.read_excel(df_path)
    df = dataframe_cleanup(df, ["Species"])

    # get the compartments and map it to the corresponding column in the df
    init_count_columns = df.filter(like='init count').columns
    compartments_to_init_map = {col.split('init count')[0].strip(): col for col in init_count_columns}

    # get the available compartments from the simulation instance
    compartments = []
    for key, value in sim_manager.simulation._children.items():  # extract the compartment types from the geometry
        if isinstance(value, stgeom._TetCompartment) or isinstance(value, stgeom._TetPatch):
            compartments.append(value)

    # set the initial values for each species in each compartment
    for compartment in compartments_to_init_map.keys():
        for s_idx, species in enumerate(sim_manager.species_names):
            # get the desired initial_value from the df
            initial_value = df.loc[df.Species == species][compartments_to_init_map[compartment]].values
            initial_value = 0 if np.isnan(initial_value) else initial_value
            try:
                # set the initial value. If this species does not exist in the compartment, ignore the error and continue
                getattr(getattr(sim_manager.simulation, compartment), species).Count = initial_value * factor
                value = getattr(getattr(sim_manager.simulation, compartment), species).Count
                print(compartment, species, value)
            except steps.API_2.sim.SolverCallError:
                pass

    
# Helperfunction to set project root and add it to path

def get_repo_path():
    """
    Retrieves the repository path based on the current user and hostname.

    This function reads the configuration from a `config.json` file located in the
    `steps_cell_signaling` directory. It dynamically determines the repository path
    based on the current system's username and hostname.

    Returns:
        str: The repository path for the current user and hostname.

    Raises:
        FileNotFoundError: If the `config.json` file is not found.
        ValueError: If no repository path is found for the current user and hostname.
    """
    # Path to the configuration file
    config_file = os.path.join(os.path.dirname(__file__), '..','..', 'config.json')

    # Load the configuration file
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file '{config_file}' not found. Please add the config.json as described i the README.")
    except json.JSONDecodeError:
        raise ValueError(f"Error parsing the configuration file '{config_file}'.")

    # Get the current hostname and username
    hostname = socket.gethostname()
    user = getpass.getuser()

    # Retrieve the repository path from the configuration
    user_config = config.get("users", {}).get(user, {})
    repo_path = user_config.get(hostname, None)

    if repo_path:
        sys.path.append(repo_path)  # Add the repository path to the Python path
        return repo_path
    else:
        raise ValueError(f"No repository path found for user '{user}' and hostname '{hostname}'.")