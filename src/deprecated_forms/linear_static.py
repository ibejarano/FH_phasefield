from dolfin import sym, grad, tr, Identity, inner, dx, TrialFunction, TestFunction

def epsilon(u):
    return sym(grad(u))

def sigma(u: TrialFunction, _lambda: float, _mu: float):
    return _lambda*tr(epsilon(u))*Identity(2) + 2*_mu*epsilon(u)

def elastic_energy_funcional(u: TrialFunction, v: TestFunction, _lambda: float, _mu: float):
    return inner(sigma(u, _lambda, _mu), epsilon(v))*dx