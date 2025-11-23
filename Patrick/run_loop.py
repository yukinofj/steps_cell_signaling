import sys
from src.SimManager import SimManager
import logging
from parameters import p
from Patrick.src.Utilities import get_repo_path
import sys
sys.path.append("/home/yukinofj/code/steps_cell_signaling")


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Starting batch simulations...")

try:
    base_path = get_repo_path()

    for i in range(1, 11):  
        ellipsoidity = i / 10
        mesh_file = f"{base_path}Patrick/meshes_ellipsoidity/ellipsoidity_{ellipsoidity:.1f}.inp"
        save_file = f"{base_path}Patrick/saved_objects/ellipsoidity_runs/testrun{ellipsoidity:.1f}"

        logging.info(f"Starting simulation for ellipsoidity {ellipsoidity:.1f}")
        sm = SimManager(
            parameters=p,
            mesh_path=mesh_file,
            save_path=save_file,
            parallel=True,
            runname=f"ellipsoidity_{ellipsoidity:.1f}",
            plot_only_run=False,
            replace=True
        )

        sm.load_model(type="small", mesh_scale=1)
        sm.run(replicats=1)
        logging.info(f"Simulation for ellipsoidity {ellipsoidity:.1f} completed.\n")

except Exception as e:
    logging.error(f"An error occurred: {e}", exc_info=True)