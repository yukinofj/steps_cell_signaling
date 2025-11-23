import numpy as np
import os
import trimesh
import sys
import pymeshfix
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from Patrick.src.MeshProcessor import create_full_mesh
from Patrick.src.Utilities import get_repo_path

base_path = get_repo_path()
print("Base path:", base_path, flush=True)

ellipsoidity = np.linspace(0,1,11)
for ellip in ellipsoidity: # loops through ellipsoidity values in 10 steps
    output_file = base_path + f"Patrick/meshes_ellipsoidity/ellipsoidity_{np.round(ellip,3)}.inp"
    print("Attempting:", output_file, flush=True)
    create_full_mesh(output_file, ellip, mesh_size_min=0.166e-6*4, mesh_size_max=0.4e-6*4) # for each ellip, creates full mesh
    print("Finished:", output_file, flush=True)

# mesh transition from sphere to ellipsoid by creating parametric mesh, where radii are scaled accordingly