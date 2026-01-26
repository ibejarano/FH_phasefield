import logging
from typing import Tuple, Dict
from dolfin import ( Mesh, Constant, inner, grad, 
    project, Function, conditional, gt, assemble, dx,
    LinearVariationalProblem, LinearVariationalSolver
)
from src.core.physics import epsilon, sigma, psi_positive
from scipy.optimize import root_scalar
from mpi4py import MPI

from src.fields.displacement import DisplacementField
from src.fields.phase import PhaseField
from src.fields.history import HistoryField
from src.fields.stress import StressField
from .config import MaterialProperties, SimulationConfig

logger = logging.getLogger(__name__)

class HydraulicFractureModel:
    def __init__(self, mesh: Mesh, material: MaterialProperties, config: SimulationConfig):
        self.mesh = mesh
        self.material = material
        self.config = config
        
        # Geometría / Parámetros de malla
        self.h_elem = mesh.hmin()
        # Physical length scale (regularization)
        self.l_c = material.l_c
        
        # Inicialización de campos (Fields)
        self.history = HistoryField(mesh, material.lmbda, material.mu)
        self.displacement = DisplacementField(mesh)
        self.phase = PhaseField(mesh)
        self.stress = StressField(mesh, material.lmbda, material.mu)
        
        # Parámetro de Presión (FEniCS Constant para la forma variacional)
        self.pressure_param = Constant(config.p_init)
        
    # --- Funciones de Física (Nomenclatura Académica) ---

    def _strain(self, u):
        return epsilon(u)

    def _stress_undegraded(self, u):
        return sigma(u, self.material.lmbda, self.material.mu)

    def _strain_energy_density_positive(self, u):
        return psi_positive(u, self.material.lmbda, self.material.mu)

    def _compute_fracture_volume(self, phi, u):
        """
        Calcula el volumen inducido por la apertura de la fractura.
        Matemáticamente: integral( -u . grad(phi) ) dx
        Representa el espacio 'hueco' creado por el desplazamiento u en la zona dañada phi.
        """
        return assemble(inner(grad(phi), -u) * dx)

    # --- Inicialización del Problema ---

    def initialize_problem(self, bcs_u, bcs_phi):
        """
        Define las ecuaciones (Formas Variacionales) y configura los Solvers de FEniCS.
        """
        # Funciones de Prueba (Trial) y Test
        p = self.phase.get_trialfunction()
        q = self.phase.get_testfunction()
        u = self.displacement.get_trialfunction()
        v = self.displacement.get_testfunction()
        
        p_old = self.phase.get_old()
        H = self.history.get()
        
        # 1. Forma Variacional del Desplazamiento (Equilibrio Mecánico)
        # Incluye la degradación (1-phi)^2 y el término de presión del fluido.
        k_res = 1.0e-6
        self.displacement_form = ((1 - p_old)**2 + k_res) * inner(self._strain(v), self._stress_undegraded(u)) * dx \
                                 + self.pressure_param * inner(v, grad(p_old)) * dx
        
        # 2. Forma Variacional del Campo de Fase (Evolución de Daño)
        # Basado en la regularización AT1/AT2 del funcional de Griffith.
        gc = self.material.Gc
        lc = self.l_c
        self.phase_field_form = (gc * lc * inner(grad(p), grad(q)) \
                                 + ((gc / lc) + 2.0 * H) * inner(p, q) \
                                 - 2.0 * H * q) * dx

        # Configuración de Solvers en los campos respectivos
        self.displacement.setup_solver(self.displacement_form, bcs_u)
        self.phase.setup_solver(self.phase_field_form, bcs_phi)
        
        # Resolución inicial para establecer el estado de fractura inicial
        self.phase.solve()
        
    # --- Orquestación de la Física ---

    def solve_time_step(self, dV: float) -> Tuple[float, float]:
        """
        Resuelve un paso de tiempo usando el algoritmo Staggered (Alternado).
        Ajusta la presión para satisfacer el incremento de volumen dV.
        """
        phi_tol = self.config.tol_phi
        err_phi = 1.0
        outer_ite = 0
        
        # Volumen acumulado previo
        v0 = self._compute_fracture_volume(self.phase.get_old(), self.displacement.get())
        vol_target = v0 + dV
        
        # Snapshot de la historia al inicio del paso (H_n)
        H_n = Function(self.history.V)
        H_n.assign(self.history.get())
        
        while err_phi > phi_tol:
            outer_ite += 1
            
            # A. Ajustar presión para cumplir con el volumen objetivo
            ite_p, pn = self._adjust_pressure(vol_target)
            if ite_p < 0:
                raise RuntimeError("El solver de presión no convergió.")
            
            # B. Actualizar Historia con el nuevo desplazamiento (si la fractura crece)
            self.history.update(self.displacement.get(), old_field=H_n)
            
            # C. Resolver campo de fase con la nueva deformación y nueva historia
            err_phi = self.phase.solve()
            
            if outer_ite > self.config.max_staggered_iter:
                raise RuntimeError(f"Bucle Staggered falló (err_phi={err_phi:.2e})")

        vol_final = self._compute_fracture_volume(self.phase.get(), self.displacement.get())
        return float(self.pressure_param), vol_final
        
    def commit_state(self):
        """Officializes the current state as the 'old' state for the next step."""
        self.displacement.update()
        self.phase.update()
        
        # History update is already cumulative, but ensuring consistency (final check)
        self.history.update(self.displacement.get())
        # The line below was redundant/incorrect manual updated
        # psi_pos = self._strain_energy_density_positive(self.displacement.get())
        # new_H = conditional(gt(psi_pos, self.history.get()), psi_pos, self.history.get())
        # self.history.field.assign(project(new_H, self.history.V))
        
        self.stress.update(self.phase.get(), self.displacement.get())

    def _adjust_pressure(self, vol_target: float) -> Tuple[int, float]:
        """
        Busca mediante el método de la secante la presión que genera el volumen vol_target.
        """
        def objective(p_val):
            self.pressure_param.assign(p_val)
            self.displacement.solve()
            u_new = self.displacement.get()
            return self._compute_fracture_volume(self.phase.get_old(), u_new) - vol_target
        
        p_guess = float(self.pressure_param)
        try:
            # Usamos secante porque no tenemos el Jacobiano explícito de la relación P-V
            res = root_scalar(
                objective, 
                x0=p_guess, 
                x1=p_guess * 1.05 if p_guess != 0 else 100.0, 
                method='secant', 
                xtol=self.config.tol_p,
                rtol=self.config.tol_p # Using same tolerance for relative check
            )
            return res.iterations, res.root
        except Exception as e:
            logger.error(f"Error en ajuste de presión: {e}")
            return -1, p_guess

    def compute_openings(self) -> Tuple[float, float]:
        """
        Calcula la apertura de la fractura (aperture) en el centro de la grieta.
        """
        from src.utils import compute_opening_overtime
        if MPI.COMM_WORLD.size == 1:
            return compute_opening_overtime(
                self.displacement.get(), 
                self.phase.get(), 
                self.h_elem / 20
            )
        return 0.0, 0.0

    def get_current_state(self) -> Dict:
        """Devuelve el estado actual de todos los campos físicos."""
        return {
            "u": self.displacement.get(),
            "phi": self.phase.get(),
            "stress": self.stress.get(),
            "history": self.history.get()
        }