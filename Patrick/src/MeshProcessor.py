import numpy as np
import trimesh
import pymeshfix
import gmsh
import sys


def fix_surface_holes(import_file, output_file):
    '''
    Fixes potential holes in the surface of a mesh.
    This was written as the meshes provided by the collaboration partners typically contain cut-off areas which result in
    holes in the surface, but meshing programs typically require watertight surface meshes for 3D meshing.
    '''

    # # Load the STL file using Trimesh, e.g.
    # import_file = "/home/pb/Downloads/bio_data_meshes/0022_0001_accelerator_20210315_bakal01_erk_main_21-03-15_12-37-27.off"
    # output_file = "/home/pb/Downloads/bio_data_meshes/fixed_test.stl"
    mesh = trimesh.load(import_file)

    # Ensure the mesh is valid and not empty
    if mesh.is_empty or mesh.faces.shape[0] == 0:
        raise ValueError("The mesh is empty or has no faces!")

    # Use pymeshfix with the vertex and face arrays
    fixer = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
    fixer.repair()
    fixer.save(output_file)



def create_ellipsoid_surface(center_x, center_y, center_z, r_a, r_b, r_c, mesh_size=0.05):
    """
    Creates an ellipsoid surface in Gmsh using the Python API.
    # creates watertight 2d surfaces (like membranes); volumes are between surfaces
    Args:
        center_x (float): X-coordinate of the ellipsoid center.
        center_y (float): Y-coordinate of the ellipsoid center.
        center_z (float): Z-coordinate of the ellipsoid center.
        r_a (float): Radius along the X-axis.
        r_b (float): Radius along the Y-axis.
        r_c (float): Radius along the Z-axis.
        mesh_size (float): Mesh size for the points.

    Returns:
        int: The tag of the created surface loop.
    """

    # Define the points
    p_center = gmsh.model.geo.addPoint(center_x, center_y, center_z, mesh_size)
    p_ax = gmsh.model.geo.addPoint(center_x + r_a, center_y, center_z, mesh_size)
    p_ay = gmsh.model.geo.addPoint(center_x, center_y + r_b, center_z, mesh_size)
    p_az = gmsh.model.geo.addPoint(center_x, center_y, center_z + r_c, mesh_size)
    p_nx = gmsh.model.geo.addPoint(center_x - r_a, center_y, center_z, mesh_size)
    p_ny = gmsh.model.geo.addPoint(center_x, center_y - r_b, center_z, mesh_size)
    p_nz = gmsh.model.geo.addPoint(center_x, center_y, center_z - r_c, mesh_size)

    # Define the ellipses create Arcs to connect the points
    l1 = gmsh.model.geo.addEllipseArc(p_ax, p_center, p_nz, p_nz)
    l2 = gmsh.model.geo.addEllipseArc(p_nz, p_center, p_nx, p_nx)
    l3 = gmsh.model.geo.addEllipseArc(p_nx, p_center, p_az, p_az)
    l4 = gmsh.model.geo.addEllipseArc(p_az, p_center, p_ax, p_ax)
    l5 = gmsh.model.geo.addEllipseArc(p_ax, p_center, p_ay, p_ay)
    l6 = gmsh.model.geo.addEllipseArc(p_ay, p_center, p_nx, p_nx)
    l7 = gmsh.model.geo.addEllipseArc(p_nx, p_center, p_ny, p_ny)
    l8 = gmsh.model.geo.addEllipseArc(p_ny, p_center, p_ax, p_ax)
    l9 = gmsh.model.geo.addEllipseArc(p_nz, p_center, p_ay, p_ay)
    l10 = gmsh.model.geo.addEllipseArc(p_ay, p_center, p_az, p_az)
    l11 = gmsh.model.geo.addEllipseArc(p_az, p_center, p_ny, p_ny)
    l12 = gmsh.model.geo.addEllipseArc(p_ny, p_center, p_nz, p_nz)

    # Define the curve loops to connect the Arcs to create watertight boundaries
    ll1 = gmsh.model.geo.addCurveLoop([l5, l10, l4])
    ll2 = gmsh.model.geo.addCurveLoop([l9, -l5, l1])
    ll3 = gmsh.model.geo.addCurveLoop([-l10, l6, l3])
    ll4 = gmsh.model.geo.addCurveLoop([-l6, -l9, l2])
    ll5 = gmsh.model.geo.addCurveLoop([l8, -l4, l11])
    ll6 = gmsh.model.geo.addCurveLoop([l12, -l8, -l1])
    ll7 = gmsh.model.geo.addCurveLoop([-l11, -l3, l7])
    ll8 = gmsh.model.geo.addCurveLoop([-l2, -l7, -l12])
    #
    # Create surfaces from the curve loops to fill watertight boundaries
    s1 = gmsh.model.geo.addSurfaceFilling([ll1])
    s2 = gmsh.model.geo.addSurfaceFilling([ll2])
    s3 = gmsh.model.geo.addSurfaceFilling([ll3])
    s4 = gmsh.model.geo.addSurfaceFilling([ll4])
    s5 = gmsh.model.geo.addSurfaceFilling([ll5])
    s6 = gmsh.model.geo.addSurfaceFilling([ll6])
    s7 = gmsh.model.geo.addSurfaceFilling([ll7])
    s8 = gmsh.model.geo.addSurfaceFilling([ll8])

    # Define the surface loop add watertight surfaces to create the watertight surfacemesh
    sl1 = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7, s8])
    # Kernel update to be available to all functions
    gmsh.model.geo.synchronize()
    return sl1


def generate_ellipsoid_radii(ellipsoidity, volume=1):
    """
    Generate an ellipsoid, increasing ellipsoidity transitions from a sphere ( = 0) to an ellipsoid ( > 0)
    while maintaining a specified constant volume.

    Parameters:
    ellipsoidity : float
        Controls how elongated the ellipsoid becomes (higher values mean more elongation).
    volume : float, optional
        The desired volume of the ellipsoid (default is 1).

    Returns:
    List of tuples (r_a, r_b, r_c) representing ellipsoid radii.
    """

    # Initial radius for a sphere with given volume
    r_0 = (3 * volume / (4 * np.pi)) ** (1 / 3)

    factor = 1 + ellipsoidity

    # Define radii while keeping volume constant
    r_a = r_0 * factor
    r_b = r_0 / np.sqrt(factor)
    r_c = r_0 / np.sqrt(factor)

    computed_volume = (4 / 3) * np.pi * r_a * r_b * r_c
    # Ensure volume remains constant
    # assert np.isclose(computed_volume, volume, atol=1e-6), "Volume constraint violated!"
    return r_a, r_b, r_c

def extrude_exo_from_cytosole_create_3D_mesh(input_file, output_file, height):

    """
    Extrudes the cytosol in a 3D mesh from an input STL file using Gmsh.

    This function takes an STL file representing a biological surface mesh, typically the output from
    fix_surface_holes(), then extrudes an extracellular layer, and generates a 3D volumetric mesh that includes the
    nucleus. The output must be saved in Gmsh's `.msh` format.

    Parameters:
    -----------
    input_file : str
        Path to the input STL file (typically the output from a surface hole-fixing process).
    output_file : str
        Path to save the output 3D mesh in `.msh` format.
    height : float
        The extrusion height (negative value) to create an extracellular layer.


    Notes:
    ------
    - The function assumes the STL file represents a closed surface.
    - The output mesh is stored in Gmsh 2.2 ASCII format for compatibility with STEPS.

    Example Usage:
    --------------
    ```python
    extrude_cytosol_include_nucleus("input.stl", "output.msh", height=5.0)
    ```
    """

    assert output_file[-4:] == ".msh", "We need to write to Gmsh 2.2 .msh file for STEPS to be able to read it, please choose the output_file accordingly"
    gmsh.initialize()
    gmsh.model.add("stl_mesh")

    # 1. Merge the STL file (loads it into Gmsh)
    gmsh.merge(input_file)

    # 2. Identify Surfaces from STL
    stl_surfaces = gmsh.model.getEntities(2)
    if not stl_surfaces:
        print("No surfaces found in the STL file, quitting...")
        gmsh.finalize()

    # 3. Create Surface Loop for STL (assuming it's a closed surface)
    stl_surface_tags = [s[1] for s in stl_surfaces]
    stl_surface_loop_tag = gmsh.model.geo.addSurfaceLoop(stl_surface_tags)
    stl_volume_tag = gmsh.model.geo.addVolume([stl_surface_loop_tag], 1)
    gmsh.model.geo.synchronize()

    # 4. Extrude the extracellular space from the cell surface mesh
    extruded_surface_tag, extruded_volume_tag = gmsh.model.geo.extrudeBoundaryLayer(stl_surfaces, heights=[-height], recombine=True)
    gmsh.model.geo.synchronize()

    # 5. Define mesh settings
    gmsh.option.setNumber("Mesh.Algorithm", 1)  # Change meshing algorithm if needed
    gmsh.option.setNumber("Mesh.MeshSizeMin", 3.0)  # Adjust as needed
    gmsh.option.setNumber("Mesh.MeshSizeMax", 8.0)  # Adjust as needed

    # 6. Generate the Mesh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # Save as Gmsh 2.2 ASCII format
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_file)
    gmsh.fltk.run()
    # Finalize Gmsh
    gmsh.finalize()
    print(f"3D mesh saved to {output_file}. Keep in mind that this function does not add a nucleus to the mesh!")


def create_full_mesh(
        output_file,
        ellipsoidity,
        cell_volume=1676e-18,
        nucleus_volume_ratio=0.2,
        extracellular_volume_offset=3e-6,
        mesh_algorithm=1,
        mesh_size_min=0.166e-6,
        mesh_size_max=0.4e-6,
):
    """
    Creates a 3D mesh of a biological cell, including nucleus, cytosol, and extracellular space,
    while maintaining specified volume constraints and ellipsoidity.

    Parameters:
    output_file (str): Path to save the generated mesh file.
    ellipsoidity (float): Controls the elongation of the cell structure.
    cell_volume (float, optional): Total volume of the cell in cubic meters (default: 1676e-18 m³).
    nucleus_volume_ratio (float, optional): Ratio of nucleus volume to total cell volume (default: 0.2).
    extracellular_volume_offset (float, optional): Offset added to extracellular radii (default: 3e-6 m).
    mesh_algorithm (int, optional): Gmsh meshing algorithm (default: 1).
    mesh_size_min (float, optional): Minimum mesh element size (default: 0.166e-6 m).
    mesh_size_max (float, optional): Maximum mesh element size (default: 0.4e-6 m).
    """

    gmsh.initialize()
    gmsh.model.add(f"mesh_ellipsoidity_{ellipsoidity}")
    rx, ry, rz = generate_ellipsoid_radii(ellipsoidity, cell_volume)

    # Generate extracellular space by extending cytosol radii
    extracellular_surface_tag = create_ellipsoid_surface(
        0, 0, 0,    # Center of the ellipsoid is at the origin
        rx + extracellular_volume_offset,
        ry + extracellular_volume_offset,
        rz + extracellular_volume_offset,
        0.05
    )

    # Generate cytosol ellipsoid surface
    cytosol_surface_tag = create_ellipsoid_surface(0, 0, 0, rx, ry, rz, 0.05)

    # Generate nucleus ellipsoid surface
    rx, ry, rz = generate_ellipsoid_radii(0, cell_volume * nucleus_volume_ratio)
    nucleus_surface_tag = create_ellipsoid_surface(0, 0, 0, rx, ry, rz, 0.05)

    # Define the volumes, more specifically the ellipsoid shells
    extracellular_volume_tag = gmsh.model.geo.addVolume([extracellular_surface_tag, cytosol_surface_tag])
    cytosol_volume_tag = gmsh.model.geo.addVolume([cytosol_surface_tag, nucleus_surface_tag])
    nucleus_volume_tag = gmsh.model.geo.addVolume([nucleus_surface_tag])

    gmsh.model.geo.synchronize()


    # Define mesh settings
    gmsh.option.setNumber("Mesh.Algorithm", mesh_algorithm)  # Choose meshing algorithm
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)  # Set minimum mesh size, the smaller the finer
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)  # Set maximum mesh size

    # Generate the mesh in 3D
    gmsh.model.mesh.generate(3)
    if output_file[-4:] != ".inp": # export
        print()
        print("Careful, suggested output file format compatible with STEPS and the rest of this software is .inp")
        print()
    gmsh.write(output_file)

    # Uncomment for visualization/debugging
    # gmsh.fltk.run()
    gmsh.finalize()