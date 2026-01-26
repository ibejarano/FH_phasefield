from dolfin import FunctionSpace, Function, project, conditional, gt
from src.core.physics import psi_positive as psi

class HistoryField:
    def __init__(self, mesh, _lambda, _mu):
        """
        V: espacio de funciones
        psi_func: función de energía (por ejemplo, psi_linear)
        E_expr: módulo de Young espacial (Expression)
        nu: coeficiente de Poisson
        data: diccionario de configuración
        """
        self.V = FunctionSpace(mesh, "DG", 0)
        self.field = Function(self.V)
        self._lambda = _lambda
        self._mu = _mu

    def update(self, u: Function, old_field: Function = None):
        """
        Actualiza el campo de historia con el desplazamiento actual u.
        Si old_field es provisto, H = max(old_field, psi(u)).
        De lo contrario, H = max(H_current, psi(u)).
        """
        psi_val = psi(u, self._lambda, self._mu)
        
        # Base field for comparison (H_n or current H)
        base_field = old_field if old_field is not None else self.field
        
        # Actualiza el campo de historia: H = max(base, psi)
        new_H = conditional(gt(psi_val, base_field), psi_val, base_field)
        self.field.assign(project(new_H, self.field.function_space()))

    def get(self):
        """
        Devuelve el campo de historia actual (Function).
        """
        return self.field