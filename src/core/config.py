from dataclasses import dataclass, field
from typing import Optional

@dataclass
class MaterialProperties:
    E: float
    nu: float
    Gc: float
    l_c: float = 0.01
    
    @property
    def mu(self):
        return self.E / (2 * (1 + self.nu))
    
    @property
    def lmbda(self):
        return self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

    def __str__(self):
        return (
            f"\nMaterial Properties:\n"
            f"  Young's Modulus (E) : {self.E:.2e} Pa\n"
            f"  Poisson's Ratio (nu): {self.nu}\n"
            f"  Fracture Energy (Gc): {self.Gc}\n"
            f"  Length Scale (l_c)  : {self.l_c}\n"
            f"  Shear Modulus (mu)  : {self.mu:.2e} Pa\n"
            f"  Lame's 1st param    : {self.lmbda:.2e} Pa"
        )

@dataclass
class SimulationConfig:
    dt: float
    t_max: float
    Q0: float
    p_init: float = 0.0
    # Geometric parameters (fixed to ensure physical consistency)
    l_init: float = 0.025   # Initial crack half-length
    w_init: float = 0.0015  # Initial crack width/thickness
    
    output_freq: int = 5
    store_freq: int = 5
    symmetric: bool = False
    case_dir: str = "output"
    mesh_name: str = "deep_fh.xml"
    
    adaptive_time: bool = False
    dt_min: float = 1e-6
    dt_max: float = 2e-2
    dt_growth: float = 1.1
    dt_shrink: float = 0.5
    dphi_max: float = 0.1
    
    # Tolerances & Iteration Limits
    tol_phi: float = 1e-6 
    tol_p: float = 1e-6
    max_staggered_iter: int = 50

    def __str__(self):
        return (
            f"\nSimulation Config:\n"
            f"  [Time Integration]\n"
            f"    dt          : {self.dt:.2e}\n"
            f"    t_max       : {self.t_max:.2e}\n"
            f"    Adaptive    : {self.adaptive_time}\n"
            f"    dt_growth   : {self.dt_growth}\n"
            f"    dt_shrink   : {self.dt_shrink}\n"
            f"\n"
            f"  [Injection]\n"
            f"    Q0          : {self.Q0:.2e}\n"
            f"    p_init      : {self.p_init}\n"
            f"\n"
            f"  [Geometry]\n"
            f"    Symmetric   : {self.symmetric}\n"
            f"    l_init      : {self.l_init}\n"
            f"    w_init      : {self.w_init}\n"
            f"    Case Dir    : {self.case_dir}\n"
            f"\n"
            f"  [Tolerances]\n"
            f"    tol_phi     : {self.tol_phi}\n"
            f"    tol_p       : {self.tol_p}\n"
            f"    max_iter    : {self.max_staggered_iter}"
        )
