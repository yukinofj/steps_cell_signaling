import steps.interface
import steps.model as stmodel   # chem species, reactions and diffusion
import steps.geom as stgeom     # compartments, surfaces and mesh based geometry
import steps.rng as strng
import steps.sim as stsim       # stochastic solver: TetOpSplit
import steps.saving as stsave   # how and where sim results stored
from src.Utilities import molar_to_molecules, nostdout
import numpy as np
import os


def initialize_ellipsoid_mesh(mesh_path, scale, nucleus_volume, cytosol_volume, extracellular_volume, cell_surface_system):
    # Load mesh and compartments
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale) # loads tetrahedral mesh

    with mesh:

        # LISTEN

        # Extracellular space
        exo_tets = stgeom.TetList(mesh.tetGroups["Volume1"])

        # Zelle
        cytosol_tets = stgeom.TetList(mesh.tetGroups["Volume2"])

        # Zellkern
        nuc_tets = stgeom.TetList(mesh.tetGroups["Volume3"])

        # # TETS im Cyt, die an die Membran grenzen
        # mem_tris = cytosol_tets.surface
        # mem_tet = stgeom.TetList()
        # for tri in mem_tris:
        #     for tet in tri.tetNeighbs:
        #         mem_tet.append(tet)
        # mem_tet = stgeom.TetList([tet for tet in mem_tet if tet in cytosol_tets])
        #
        # # Hälfte Cyto
        # half_cyt = stgeom.TetList(tet for tet in cytosol_tets if tet.center.y > 0)
        # # TRIS vom Durchschnitt
        # slice_cyt = stgeom.TriList(half_cyt.surface & (cytosol_tets - half_cyt).surface)
        # # die Tets am Querschnitt
        # slice_tet_cyt = stgeom.TetList()
        # triInds_cyt = []
        # for tri in slice_cyt:
        #     for tet in tri.tetNeighbs:
        #         slice_tet_cyt.append(tet)
        #         triInds_cyt.append(tri.idx)
        #
        # # Häfte Nuc
        # half_nuc = stgeom.TetList(tet for tet in nuc_tets if tet.center.y > 0)
        # # TRIS vom Durchscnitt
        # slice_nuc = stgeom.TriList(half_nuc.surface & (nuc_tets - half_nuc).surface)
        # # die Tets am Querschnitt
        # slice_tet_nuc = stgeom.TetList()
        # triInds_nuc = []
        # for tri in slice_nuc:
        #     for tet in tri.tetNeighbs:
        #         slice_tet_nuc.append(tet)
        #         triInds_nuc.append(tri.idx)


        # COMPARTMENTS
        # Zellkern
        nuc = stgeom.Compartment(nuc_tets, nucleus_volume, name="nuc")

        # Cytoplasma
        cyt = stgeom.Compartment(cytosol_tets, cytosol_volume, name="cyt")

        # Zelläußeres
        exo = stgeom.Compartment(exo_tets, extracellular_volume, name="exo")

        # Zellmembran: between cytosol and extracellular space, membrane reactions occur here
        cell_surface = stgeom.Patch(cyt.surface & exo.surface, cyt, exo, cell_surface_system, name="cell_surface")

        # DIFFUSIONS BARRIERE, why is this necessary anywhere? We have discrete volumes anyways, so I dont think this makes sense
        # Zellkernmembran
        nuc_mem = stgeom.DiffBoundary(nuc.surface, name="nuc_mem") # barrier for nuc
    return mesh, exo_tets, cytosol_tets, nuc_tets


def create_model(p, species_names, mesh_path, mesh_scale, plot_only_run):
    """
    Creates a STEPS simulation model based on parameters, species, and mesh geometry.

    Parameters:
    -----------
    p : dict
        Simulation parameters including diffusion constants (`p["DC"]`), reaction rates,
        and time step (`p["time step"]`).

    species_names : list of str
        Names of chemical species to be defined in the model.

    mesh_path : str
        Path to the Abaqus mesh file defining system geometry, scaled to micrometers.

    Returns:
    --------
    simulation : steps.sim.Simulation
        Configured STEPS simulation object.

    rs : steps.saving.ResultSelector
        Object for selecting and storing simulation results.

    Raises:
    -------
    AssertionError:
        If the provided `mesh_path` is invalid.

    Notes:
    ------
    - Defines species diffusion, surface and volumetric reactions (e.g., EGF-EGFR binding).
    - Partitions the mesh into compartments: nucleus, cytoplasm, extracellular space,
      and membrane patch.
    - Sets up a stochastic simulation engine with result selectors for species counts.
    """
    mdl = stmodel.Model() # global model
    r = stmodel.ReactionManager()
    species_dict = {}

    # Create volume and surface systems
    with mdl: # reaction diffusion environments
        cytosol_volume = stmodel.VolumeSystem(name="cytosol_volume")
        nucleus_volume = stmodel.VolumeSystem(name="nucleus_volume")
        extracellular_volume = stmodel.VolumeSystem(name="extracellular_volume")
        cell_surface_system = stmodel.SurfaceSystem(name = "cell_surface")

        # Create a dictionary to hold the created species
        for sp_name in species_names:
            species_dict[sp_name] = stmodel.Species(name=sp_name)

        # Cytoplasma/Cytosol
        with cytosol_volume:
            # ERK Deaktivierung (dephosphorylation)
            species_dict["ERKp"] + species_dict["P3"] > r[6] > species_dict["ERK"] + species_dict["P3"]
            r[6].K = p["k[666]"]
            #Cytoplasm diffusion
            stmodel.Diffusion(species_dict["ERK"], p["DC"])
            stmodel.Diffusion(species_dict["ERKp"], p["DC"])
            stmodel.Diffusion(species_dict["P3"], p["DC"] * 2) # faster
            stmodel.Diffusion(species_dict["GAP"], p["DC"]/4) # slower

        # Extracellular volume
        with extracellular_volume: # free diffusion of EGF
            stmodel.Diffusion(species_dict["EGF"], p["DC"]/10)

        # Nucleus volume
        with nucleus_volume:
            species_dict["ERKp"] > r[7] > species_dict["ERKpp"] # phosphorylation
            r[7].K = 1e8
            stmodel.Diffusion(species_dict["ERKpp"], p["DC"]) # ERKpp moves within nucleus

        # Surface system (cell membrane)
        with cell_surface_system: # EGF-EGFR signaling cascade
            species_dict["EGFR"].s + species_dict["EGF"].o < r[1] > species_dict["EGF_EGFR"].s                                              # EFG binds EGFR
            species_dict["EGF_EGFR"].s + species_dict["EGF_EGFR"].s < r[2] > species_dict["EGF_EGFR2"].s                                    # receptor dimerizes
            species_dict["EGF_EGFR2"].s < r[3] > species_dict["EGF_EGFRp2"].s                                                               # phosphorylation
            species_dict["EGF_EGFRp2"].s + species_dict["GAP"].i < r[4] > species_dict["EGF_EGFRp2_GAP"].s                                  # binds GAP
            species_dict["EGF_EGFRp2_GAP"].s + species_dict["ERK"].i < r[5] > species_dict["EGF_EGFRp2_GAP"].s + species_dict["ERKp"].i     # activates ERK to ERKp

            # reaction rates
            r[1].K = 3e7, 38e-4   # 1/Ms
            r[2].K = 1e7 * 100, 0.1   # 1/Ms
            r[3].K = 1 * 100, 0.01  # 1/s
            r[4].K = 1e6 * 100 , 0.2   # 1/Ms 1e6, 0.2
            r[5].K = p["k[0]"], 0.1 #1e8 * c1, 0.1 * c1

            stmodel.Diffusion(species_dict["EGFR"], p["DC"]/10) # slow diffusion
            stmodel.Diffusion(species_dict["EGF_EGFR"], p["DC"]/20)
            stmodel.Diffusion(species_dict["EGF_EGFR2"], p["DC"]/40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2"], p["DC"]/40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP"], p["DC"]/50)

    # Load mesh and compartments
    mesh, exo_tets, cytosol_tets, nuc_tets = initialize_ellipsoid_mesh(mesh_path,
                                                                    scale=mesh_scale,
                                                                    nucleus_volume=nucleus_volume,
                                                                    cytosol_volume=cytosol_volume,
                                                                    extracellular_volume=extracellular_volume,
                                                                    cell_surface_system=cell_surface_system)
    system_volume = mesh.Vol

    # Initialize RNG and Simulation
    rng = strng.RNG("mt19937", 512, 2903) # (type of RNG (mt: mersenne twister), number of pregenerated random numbers, seed value to initialize the RNG)
    # with nostdout(): #doesnt work, crashes the freaking sim...
    partition = stgeom.LinearMeshPartition(mesh, 1, 1, stsim.MPI.nhosts)
    simulation = stsim.Simulation("TetOpSplit", mdl, mesh, rng, False, partition)


    # Define which results to save and how they are processed
    # Not all species are present in all the compartments, be aware of this.
    rs = None
    if plot_only_run == False:
        rs = stsave.ResultSelector(simulation)
        result_selectors = {
            "EGF": rs.SUM(rs.TETS(exo_tets).EGF.Count),
            "EGF_EGFR": rs.SUM(rs.TRIS(cytosol_tets.surface).EGF_EGFR.Count),
            "EGF_EGFR2": rs.SUM(rs.TRIS(cytosol_tets.surface).EGF_EGFR2.Count),
            "EGF_EGFRp2": rs.SUM(rs.TRIS(cytosol_tets.surface).EGF_EGFRp2.Count),
            "EGF_EGFRp2_GAP": rs.SUM(rs.TRIS(cytosol_tets.surface).EGF_EGFRp2_GAP.Count),
            "ERK_cyto": rs.SUM(rs.TETS(cytosol_tets).ERK.Count),
            "ERKp_cyto": rs.SUM(rs.TETS(cytosol_tets).ERKp.Count),
            "ERKp_nuc": rs.SUM(rs.TETS(nuc_tets).ERKp.Count),
            "ERKpp": rs.SUM(rs.TETS(nuc_tets).ERKpp.Count),
            # "Concentrations": rs.TETS().LIST(species_dict["EGF_EGFR"],
            #                                  species_dict["ERK"],
            #                                  species_dict["ERKpp"]).Conc,
        }
        # Schedule saving
        for key, sel in result_selectors.items():
            simulation.toSave(sel, dt=p["time step"])

    return simulation, rs, mesh
