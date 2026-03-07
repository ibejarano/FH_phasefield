import argparse
from src.core.config import MaterialProperties, SimulationConfig
from src.core.runner import run_simulation
from src.utils import setup_logging

setup_logging()

parser = argparse.ArgumentParser(description="Radial Transition Hydraulic Fracture Simulation (axisymmetric)")
parser.add_argument("case_name", nargs="?", default="output_transition", help="Output directory name")
parser.add_argument("--mesh", type=str, default="geo_files/transition_axi_unstruc.geo", help="Path to .geo mesh file")
parser.add_argument("--debug", action="store_true", help="Enable debug output for convergence monitoring")
args = parser.parse_args()

mat_props = MaterialProperties(E=2e8, nu=0.3, Gc=1)

sim_config = SimulationConfig(
    dt=1e-7,
    t_max=0.2,
    Q0=0.5e-3,
    p_init=1000,
    l_init=0.01,
    store_freq=20,
    output_freq=20,
    case_dir=args.case_name,
    symmetric=False,
    axisymmetric=True,
    tol_p=1e-8,
    tol_phi=1e-3,
    adaptive_time=True,
    dt_growth=1.03,
    dt_shrink=0.9,
    dt_min=1e-9,
    dt_max=1e-1,
    max_staggered_iter=15,
    debug=args.debug,
    sigma_h=1e2,
    phi_anisotropy=1.0,
    phi_band_width=0.0021,
)

run_simulation(mat_props, sim_config, args.mesh, upper_face_free=True)
