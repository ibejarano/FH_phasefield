#!/usr/bin/env python
"""
Generic runner for hydraulic fracture simulations.

Usage:
    python run_case.py configs/kgd_2d.json
    python run_case.py configs/kgd_2d.json --case-name output_kgd_1 --debug
    python run_case.py configs/radial_2d.json --mesh geo_files/other.geo
    python run_case.py configs/kgd_2d.json --dry-run
"""
import argparse
import json
import sys

from src.core.config import MaterialProperties, SimulationConfig
from src.core.runner import run_simulation
from src.utils import setup_logging


def load_config(json_path: str) -> dict:
    """Load and return a JSON configuration file."""
    with open(json_path, "r") as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Run a hydraulic fracture simulation from a JSON config file."
    )
    parser.add_argument("config", type=str, help="Path to the JSON config file")
    parser.add_argument("--case-name", type=str, default=None,
                        help="Override output directory name (default: from JSON)")
    parser.add_argument("--mesh", type=str, default=None,
                        help="Override mesh .geo file path (default: from JSON)")
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug output for convergence monitoring")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the configuration and exit without running")
    args = parser.parse_args()

    # Load JSON config
    cfg = load_config(args.config)

    # Build MaterialProperties
    mat_props = MaterialProperties.from_dict(cfg.get("material", {}))

    # Build SimulationConfig with optional CLI overrides
    sim_dict = cfg.get("simulation", {})
    if args.case_name is not None:
        sim_dict["case_dir"] = args.case_name
    if args.debug:
        sim_dict["debug"] = True
    sim_config = SimulationConfig.from_dict(sim_dict)

    # Mesh path: CLI override > JSON > default
    mesh = args.mesh or cfg.get("mesh", "geo_files/sym_deep.geo")

    # Runner kwargs (upper_face_free, auto_l_init_factor)
    runner_opts = cfg.get("runner", {})
    upper_face_free = runner_opts.get("upper_face_free", False)
    auto_l_init_factor = runner_opts.get("auto_l_init_factor", 0.0)

    if args.dry_run:
        print(f"Config file : {args.config}")
        print(f"Description : {cfg.get('description', '(none)')}")
        print(f"Mesh        : {mesh}")
        print(mat_props)
        print(sim_config)
        print(f"\n  [Runner Options]")
        print(f"    upper_face_free    : {upper_face_free}")
        print(f"    auto_l_init_factor : {auto_l_init_factor}")
        sys.exit(0)

    setup_logging()
    run_simulation(mat_props, sim_config, mesh,
                   upper_face_free=upper_face_free,
                   auto_l_init_factor=auto_l_init_factor)


if __name__ == "__main__":
    main()
