from dolfin import Function, FunctionSpace, TrialFunction, TestFunction, errornorm
import logging
import time

logger = logging.getLogger(__name__)


class PhaseField:
    def __init__(self, mesh):
        self.mesh = mesh
        self.V = FunctionSpace(mesh, "CG", 1)
        self.new = Function(self.V, name="phi")
        self.old = Function(self.V)

    def update(self):
        self.old.assign(self.new)

    def get_trialfunction(self):
        return TrialFunction(self.V)

    def get_testfunction(self):
        return TestFunction(self.V)

    def get(self):
        return self.new
    
    def get_old(self):
        return self.old
    
    def get_error(self):
        """
        Returns the relative L2 error norm: ||new - old|| / (||new|| + epsilon)
        """
        from dolfin import norm
        diff = errornorm(self.new, self.old, norm_type='l2', mesh=self.mesh)
        norm_new = norm(self.new, norm_type='l2', mesh=self.mesh)
        return diff / (norm_new + 1.0e-12)
    
    def setup_solver(self, E_phi, bc_phi):
        """
        Configura el Solver para el campo de fase.
        """
        from dolfin import LinearVariationalProblem, LinearVariationalSolver, lhs, rhs
        p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), self.new, bc_phi)
        solver_phi = LinearVariationalSolver(p_phi)
        #solver_phi.parameters["linear_solver"] = "gmres"
        #solver_phi.parameters["preconditioner"] = "ilu"
        self.solver = solver_phi

    def solve(self):
        """
        Resuelve el problema de campo de fase.
        """
        logger.debug("Solving phase field problem...")
        start_time = time.time()
        self.solver.solve()
        error = self.get_error()
        self.update()
        elapsed_time = time.time() - start_time
        logger.debug(f"Phase field solved in {elapsed_time:.4f} seconds with error: {error:.6f}")
        return error