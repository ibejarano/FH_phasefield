import os
import logging
from mpi4py import MPI
from dolfin import Mesh, set_log_level, LogLevel

from .config import MaterialProperties, SimulationConfig
from .model import HydraulicFractureModel
from .simulator import HydraulicSimulator
from .boundaries import setup_boundary_conditions
from .meshing import generate_mesh

logger = logging.getLogger(__name__)


def run_simulation(
    mat_props: MaterialProperties,
    sim_config: SimulationConfig,
    geo_file: str,
    upper_face_free: bool = False,
    auto_l_init_factor: float = 0.0,
):
    """
    Entry point compartido: genera la malla, inicializa el modelo y ejecuta la simulación.

    Args:
        mat_props: Propiedades del material.
        sim_config: Configuración de la simulación (case_dir, symmetric, axisymmetric, etc.).
        geo_file: Ruta al archivo .geo de la malla.
        upper_face_free: Si True, la cara superior es libre (tracción libre). Default False.
        auto_l_init_factor: Si > 0, sobreescribe sim_config.l_init = mesh.hmin() * factor.
    """
    comm = MPI.COMM_WORLD
    set_log_level(LogLevel.ERROR)

    if comm.rank == 0:
        os.makedirs(sim_config.case_dir, exist_ok=True)
        _setup_file_logging(sim_config, geo_file, mat_props)
        generate_mesh(geo_file, sim_config.case_dir, "mesh")

    comm.Barrier()

    mesh = Mesh(comm, os.path.join(sim_config.case_dir, "mesh.xml"))
    sim_config.w_init = mesh.hmin()
    if auto_l_init_factor > 0.0:
        sim_config.l_init = mesh.hmin() * auto_l_init_factor

    model = HydraulicFractureModel(mesh, mat_props, sim_config)

    bcs_u, bcs_phi = setup_boundary_conditions(
        model.phase,
        model.displacement,
        l_init=sim_config.l_init,
        h_elem=sim_config.w_init,
        symmetric=sim_config.symmetric,
        axisymmetric=sim_config.axisymmetric,
        upper_face_free=upper_face_free,
        phi_band_width=getattr(sim_config, "phi_band_width", 0.0),
    )

    model.initialize_problem(bcs_u, bcs_phi)
    HydraulicSimulator(model, sim_config).run()


def _setup_file_logging(sim_config: SimulationConfig, geo_file: str, mat_props: MaterialProperties):
    log_file = os.path.join(sim_config.case_dir, "simulation.log")
    handler = logging.FileHandler(log_file)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
    logging.getLogger().addHandler(handler)

    logging.info("=" * 40)
    logging.info(f"Starting Simulation: {sim_config.case_dir}")
    logging.info("=" * 40)
    logging.info(f"  Mesh: {geo_file}")
    logging.info(f"  Symmetric: {sim_config.symmetric} | Axisymmetric: {sim_config.axisymmetric}")
    logging.info("-" * 20)
    logging.info(f"{mat_props}")
    logging.info("-" * 20)
    logging.info(f"{sim_config}")
    logging.info("=" * 40)
