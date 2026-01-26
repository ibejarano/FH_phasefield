from dolfin import Function, VectorFunctionSpace, TrialFunction, TestFunction
import logging
import time
from mpi4py import MPI


logger = logging.getLogger(__name__)


class DisplacementField:
    def __init__(self, mesh, V=None):
        self.V = V if V is not None else VectorFunctionSpace(mesh, "CG", 1)
        self.new = Function(self.V, name="displacement")
        self.old = Function(self.V)
        self.temp = Function(self.V)

    def update(self):
        self.old.assign(self.new)

    def get_trialfunction(self):
        return TrialFunction(self.V)

    def get_testfunction(self):
        return TestFunction(self.V)

    def get(self):
        return self.new
    
    def setup_solver(self, E_du, bc_u):
        """
        Configura el Solver para el campo de desplazamiento.
        """
        from dolfin import LinearVariationalProblem, LinearVariationalSolver, lhs, rhs
        p_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), self.new, bc_u)
        solver_u = LinearVariationalSolver(p_u)
        #solver_u.parameters["linear_solver"] = "gmres"
        #solver_u.parameters["preconditioner"] = "ilu"
        self.solver = solver_u

    def solve(self):
        """
        Resuelve el problema de desplazamiento.
        """
        if logger.isEnabledFor(logging.DEBUG) and MPI.COMM_WORLD.rank == 0:
            logger.debug("Solving displacement field problem...")
            start_time = time.time()

        self.solver.solve()
        self.update()
        
        if logger.isEnabledFor(logging.DEBUG) and MPI.COMM_WORLD.rank == 0:
            elapsed_time = time.time() - start_time
            logger.debug(f"Displacement field solved in {elapsed_time:.4f} seconds")