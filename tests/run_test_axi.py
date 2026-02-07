#!/usr/bin/env python
"""
Test script for axisymmetric implementation with debug output.
Runs a short simulation (5-10 steps) to verify convergence.
"""
import sys
import os

# Add parent dir to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dolfin import Mesh, set_log_level, LogLevel
from mpi4py import MPI

from src.core.config import MaterialProperties, SimulationConfig
from src.core.model import HydraulicFractureModel
from src.core.simulator import HydraulicSimulator
from src.core.boundaries import setup_boundary_conditions
from src.core.meshing import generate_mesh
from src.utils import setup_logging

# Setup
setup_logging()
set_log_level(LogLevel.ERROR)

comm = MPI.COMM_WORLD

# Configuration - short run for testing

sim_config = SimulationConfig(
    dt=1e-7,
    t_max=5e0, 
    Q0=0.5e-3,
    p_init=1000,
    l_init=0.08,
    w_init=0.001,
    store_freq= 100,
    output_freq= 10,
    case_dir="tests/output_axi_test_2",
    axisymmetric=True,
    tol_p=1e-6,
    tol_phi=1e-3,
    adaptive_time=True,  # Fixed dt for testing
    dt_min=1e-8,
    dt_max=1e-4,
    max_staggered_iter=20,
    debug=True  # Enable debug output
)

# Mesh
if comm.rank == 0:
    os.makedirs(sim_config.case_dir, exist_ok=True)
    geo_file = "tests/cases/simple_axi.geo"
    generate_mesh(geo_file, sim_config.case_dir, "mesh")

comm.Barrier()

mesh_path = os.path.join(sim_config.case_dir, "mesh.xml")
mesh = Mesh(comm, mesh_path)

mat_props = MaterialProperties(
    E=2e8,
    nu=0.3,
    Gc=1,
    l_c=mesh.hmin()*4
)

sim_config.w_init = mesh.hmin()*1
sim_config.l_init = mesh.hmin()*10

# Model setup
model = HydraulicFractureModel(mesh, mat_props, sim_config)

bcs_u, bcs_phi = setup_boundary_conditions(
    model.phase,
    model.displacement,
    l_init=sim_config.l_init,
    h_elem=sim_config.w_init,
    axisymmetric=sim_config.axisymmetric
)

model.initialize_problem(bcs_u, bcs_phi)

# Run
print("\n" + "="*60)
print("AXISYMMETRIC TEST - Debug Mode")
print("="*60)
print(f"Mesh elements: {mesh.num_cells()}")
print(f"Time steps: {int(sim_config.t_max / sim_config.dt)}")
print("="*60 + "\n")

simulator = HydraulicSimulator(model, sim_config)
simulator.run()

print("\n" + "="*60)
print("TEST COMPLETE")
print("="*60)
