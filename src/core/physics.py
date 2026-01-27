from dolfin import sym, grad, tr, Identity, inner, dev, as_tensor

def epsilon(u):
    """Infinitesimal strain tensor."""
    return sym(grad(u))

def epsilon_axi(u, r):
    """
    Axisymmetric strain tensor for 2D domain (r, z).
    u[0] = u_r (radial displacement)
    u[1] = u_z (axial displacement)
    r    = radial coordinate (x[0])
    
    Returns a 3x3 tensor-like structure or simply the relevant components 
    packaged such that tr() and inner() work as expected for energy density.
    
    Structure:
    [ eps_rr   0      eps_rz ]
    [ 0        eps_ee 0      ]  <-- eps_theta_theta = u_r / r
    [ eps_rz   0      eps_zz ]
    """
    return as_tensor([
        [u[0].dx(0), 0, 0.5*(u[0].dx(1) + u[1].dx(0))],
        [0, u[0] / (r + 1.0e-12), 0],
        [0.5*(u[0].dx(1) + u[1].dx(0)), 0, u[1].dx(1)]
    ])
def sigma(u, lmbda, mu):
    """Linear elastic stress tensor (undegraded)."""
    eps = epsilon(u)
    return lmbda * tr(eps) * Identity(2) + 2 * mu * eps

def sigma_axi(u, r, lmbda, mu):
    """Stress tensor for axisymmetry (returns 3x3)."""
    eps = epsilon_axi(u, r)
    return lmbda * tr(eps) * Identity(3) + 2 * mu * eps
def psi_positive(u, lmbda, mu):
    """Positive part of the strain energy density (tension only)."""
    eps = epsilon(u)
    return _psi_positive_tensor(eps, lmbda, mu)

def psi_positive_axi(u, r, lmbda, mu):
    """Energy density for axisymmetry."""
    eps = epsilon_axi(u, r)
    return _psi_positive_tensor(eps, lmbda, mu)

def _psi_positive_tensor(eps, lmbda, mu):
    """Core calculation for psi positive given a strain tensor."""
    tr_eps = tr(eps)
    tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
    return 0.5 * (lmbda + mu) * tr_eps_pos**2 + mu * inner(dev(eps), dev(eps))

def elastic_energy_functional(u, v, lmbda, mu):
    """Bilinear form for linear elastic energy."""
    from dolfin import inner, dx
    return inner(sigma(u, lmbda, mu), epsilon(v)) * dx