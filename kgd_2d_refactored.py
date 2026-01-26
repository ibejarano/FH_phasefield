import sys
import os
from dolfin import Mesh, set_log_level, LogLevel
from mpi4py import MPI

# Use the refined core components
from src.core.config import MaterialProperties, SimulationConfig
from src.core.model import HydraulicFractureModel
from src.core.simulator import HydraulicSimulator
from src.core.boundaries import setup_boundary_conditions
from src.core.meshing import generate_mesh
from src.utils import setup_logging

# Setup environment
setup_logging() # Configures the root logger for all modules
set_log_level(LogLevel.ERROR)

import argparse

# 1. Configuration
parser = argparse.ArgumentParser(description="Hydraulic Fracture Simulation")
parser.add_argument("case_name", nargs="?", default="output_legible", help="Output directory name")
parser.add_argument("--symmetric", action="store_true", help="Use symmetric formulation")
parser.add_argument("--mesh", type=str, help="Path to .geo mesh file")
args = parser.parse_args()

case_name = args.case_name
is_symmetric = args.symmetric

mat_props = MaterialProperties(
    E = 2e8,
    nu = 0.3,
    Gc = 1
)

sim_config = SimulationConfig(
    dt = 2e-6,
    t_max = 1e-2    ,
    Q0 = 0.5e-3,
    p_init = 1000,
    case_dir = case_name,
    symmetric = is_symmetric,
    tol_p=1e-7,
    tol_phi=1e-3,
    adaptive_time=True,
    dt_growth=1.1,
    dt_shrink=0.5,
    dt_min=1e-6,
    dt_max=1e-2,
    max_staggered_iter=15
)

# 2. Mesh setup
comm = MPI.COMM_WORLD

# Configure File Logging & Create Directory (Rank 0 only)
if comm.rank == 0:
    import logging
    os.makedirs(sim_config.case_dir, exist_ok=True)
    
    # Add File Handler
    log_file = os.path.join(sim_config.case_dir, "simulation.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    
    # Log Parameters
    logging.info("="*40)
    logging.info(f"Starting Simulation: {case_name}")
    logging.info("="*40)
    logging.info(f"Mesh Configuration:")
    logging.info(f"  - Input Mesh: {args.mesh if args.mesh else 'Default (Symmetric/Full)'}")
    logging.info(f"  - Symmetric Mode: {is_symmetric}")
    
    logging.info("-" * 20)
    logging.info(f"{mat_props}")
    
    logging.info("-" * 20)
    logging.info(f"{sim_config}")
    logging.info("="*40)

# Only rank 0 generates the mesh
if comm.rank == 0:
    if args.mesh:
        geo_file = args.mesh
    else:
        geo_file = "deep_fracture/sym_deep.geo" if is_symmetric else "deep_fracture/deep_fh.geo"
    generate_mesh(geo_file, sim_config.case_dir, "mesh")

# Wait for mesh generation to complete
comm.Barrier()

mesh_path = os.path.join(sim_config.case_dir, "mesh.xml")
mesh = Mesh(comm, mesh_path)

# 3. Model & Boundary Conditions
model = HydraulicFractureModel(mesh, mat_props, sim_config)

# Setup BCs (using newly refactored core component) # l_init = model.l_c * 2.5
# Setup BCs (using newly refactored core component)
bcs_u, bcs_phi = setup_boundary_conditions(
    model.phase, 
    model.displacement, 
    l_init=sim_config.l_init, 
    h_elem=sim_config.w_init, # Use w_init as the effective initial width
    symmetric=is_symmetric
)

model.initialize_problem(bcs_u, bcs_phi)

# 4. Engine execution
simulator = HydraulicSimulator(model, sim_config)
simulator.run()
