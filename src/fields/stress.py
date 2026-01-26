from dolfin import Function, TensorFunctionSpace, project
import logging
from src.core.physics import sigma

logger = logging.getLogger(__name__)


class StressField:
    def __init__(self, mesh, _lambda, _mu, V=None):
        self.V = V if V is not None else TensorFunctionSpace(mesh, "DG", 0)
        self.current = Function(self.V, name="stress")
        self._lambda = _lambda
        self._mu = _mu

    def update(self, pnew: Function, unew: Function):
        stress_expr = (1-pnew)**2 * sigma(unew, self._lambda, self._mu)
        self.current.assign(project(stress_expr, self.V))

    def get(self):
        return self.current