from dolfin import sym, grad, tr, Identity, inner, dev

def epsilon(u):
    """Infinitesimal strain tensor."""
    return sym(grad(u))
def sigma(u, lmbda, mu):
    """Linear elastic stress tensor (undegraded)."""
    eps = epsilon(u)
    return lmbda * tr(eps) * Identity(2) + 2 * mu * eps
def psi_positive(u, lmbda, mu):
    """Positive part of the strain energy density (tension only)."""
    eps = epsilon(u)
    tr_eps = tr(eps)
    tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
    return 0.5 * (lmbda + mu) * tr_eps_pos**2 + mu * inner(dev(eps), dev(eps))

def elastic_energy_functional(u, v, lmbda, mu):
    """Bilinear form for linear elastic energy."""
    from dolfin import inner, dx
    return inner(sigma(u, lmbda, mu), epsilon(v)) * dx