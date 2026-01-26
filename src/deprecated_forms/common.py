from dolfin import inner, grad, dx, assemble, TrialFunction, tr, Identity, sym, dev, Function


def epsilon(u):
    return sym(grad(u))

def sigma(u: TrialFunction, _lambda: float, _mu: float):
    return _lambda*tr(epsilon(u))*Identity(2) + 2*_mu*epsilon(u)

def psi(u: Function, _lambda: float, _mu: float):
    return 0.5*(_lambda + _mu)*(0.5*(tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + _mu*inner(dev(epsilon(u)), dev(epsilon(u)))


def compute_fracture_volume(phi, u):
    vol_frac = assemble( inner(grad(phi), -u) * dx )
    return vol_frac