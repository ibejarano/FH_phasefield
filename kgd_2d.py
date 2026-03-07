import argparse
from src.core.config import MaterialProperties, SimulationConfig
from src.core.runner import run_simulation
from src.utils import setup_logging

setup_logging()

parser = argparse.ArgumentParser(description="KGD Hydraulic Fracture Simulation (plane strain, symmetric)")
parser.add_argument("case_name", nargs="?", default="output_kgd", help="Output directory name")
parser.add_argument("--mesh", type=str, default="geo_files/sym_deep.geo", help="Path to .geo mesh file")
parser.add_argument("--debug", action="store_true", help="Enable debug output for convergence monitoring")
args = parser.parse_args()

mat_props = MaterialProperties(E=2e8, nu=0.3, Gc=1)

sim_config = SimulationConfig(
    dt=1e-4,
    t_max=8,
    Q0=0.5e-3,
    p_init=1000,
    l_init=0.2,
    store_freq=20,
    output_freq=20,
    case_dir=args.case_name,
    symmetric=True,
    axisymmetric=False,
    tol_p=1e-7,
    tol_phi=5e-4,
    adaptive_time=True,
    dt_growth=1.05,
    dt_shrink=0.5,
    dt_min=1e-9,
    dt_max=1e-1,
    max_staggered_iter=25,
    debug=args.debug,
    sigma_h=1e4,
)

run_simulation(mat_props, sim_config, args.mesh, upper_face_free=True)
