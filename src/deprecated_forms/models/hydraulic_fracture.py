from dolfin import Constant, Mesh
from typing import Dict, Any, Tuple
import logging

from src.fields.displacement import DisplacementField
from src.fields.phase import PhaseField
from src.fields.history import HistoryField
from src.fields.stress import StressField
from src.variational_forms.phase_field import phase_field_problem, solve_step_staggered

logger = logging.getLogger(__name__)

class HydraulicFractureModel:
    def __init__(self, mesh: Mesh, params: Dict[str, Any]):
        self.mesh = mesh
        self.params = params
        
        # Lame parameters
        self.mu = params['E'] / (2 * (1 + params['nu']))
        self.lmbda = params['E'] * params['nu'] / ((1 + params['nu']) * (1 - 2 * params['nu']))
        self.Gc = params['Gc']
        
        # Mesh parameters
        self.h_elem = mesh.hmin()
        self.l_c = self.h_elem * params.get('relacion_hl', 4)
        
        # Fields
        self.history = HistoryField(mesh, self.lmbda, self.mu)
        self.displacement = DisplacementField(mesh)
        self.phase = PhaseField(mesh)
        self.stress = StressField(mesh, self.lmbda, self.mu)
        
        self.pressure_val = params.get('p_init', 0.0)
        self.pressure = None # Set later in setup_variational_forms
        
    def setup_problem(self, bcs_u, bcs_phi):
        """Sets up variational forms and solvers."""
        E_du, E_phi, self.pressure = phase_field_problem(
            self.phase, self.displacement, self.history,
            self.lmbda, self.mu, self.Gc,
            self.pressure_val, self.l_c
        )
        
        self.displacement.setup_solver(E_du, bcs_u)
        self.phase.setup_solver(E_phi, bcs_phi)
        
        # Initial phase solve if needed
        self.phase.solve()
        
    def solve_step(self, dV: float) -> Tuple[float, float]:
        """Solves a single time step using staggered approach."""
        pn, vol_frac = solve_step_staggered(
            self.displacement, 
            self.phase, 
            self.history, 
            self.pressure, 
            dV
        )
        
        # Update stress field for post-processing
        self.stress.update(self.phase.get(), self.displacement.get())
        
        return pn, vol_frac

    def get_fields(self):
        return {
            "u": self.displacement.get(),
            "phi": self.phase.get(),
            "stress": self.stress.get(),
            "history": self.history.get()
        }
