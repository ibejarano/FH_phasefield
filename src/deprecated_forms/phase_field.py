from dolfin import inner, grad, dx, Constant
from variational_forms.common import compute_fracture_volume, sigma, epsilon
from scipy.optimize import root_scalar

import logging

logger = logging.getLogger(__name__)

def phase_field_problem(
        phase,
        displacement,
        history, 
        lmbda: float,
        mu: float,
        Gc: float,
        p_init: float,
        l_c: float
):
    pressure = Constant(p_init)
    p = phase.get_trialfunction()
    q = phase.get_testfunction()
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    pold = phase.get_old()

    sigma_u = sigma(u, lmbda, mu)

    H = history.get()

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma_u) * dx \
           + pressure * inner(v, grad(pold)) * dx

    E_phi = (Gc * l_c * inner(grad(p), grad(q)) \
             + ((Gc / l_c) + 2.0 * H) * inner(p, q) \
             - 2.0 * H * q) * dx

    return E_du, E_phi, pressure

def compute_pressure(displacement, 
                     phase, 
                     history, 
                     pressure: float, 
                     vol_target: float, 
                     method="root_scalar"):
    
    ite_p, pn = pressure_solver(
        Vtarget=vol_target,
        phase=phase,
        displacement=displacement,
        history=history,
        pressure=pressure,
        vol_tol=1e-6,
        method=method
    )

    return ite_p, pn

def solve_step_staggered(displacement, phase, history, pressure, dV:float, phi_tol=1e-3):
    err_phi = 1.0
    outer_ite = 0
    V0 = compute_fracture_volume(phase.get_old(), displacement.get())
    vol_target = V0 + dV

    while err_phi > phi_tol:
        outer_ite += 1
        ite_p, pn = compute_pressure(displacement, phase, history, pressure, vol_target)
        if ite_p < 0:
            raise RuntimeError("Pressure adjustment failed to converge.")
        err_phi = phase.solve()
        if outer_ite > 15:
            raise RuntimeError(f"Outer staggered loop failed to converge (err_phi={err_phi:.2e})")

    displacement.update()
    phase.update()
    history.update(displacement.get())

    vol = compute_fracture_volume(phase.get(), displacement.get())

    return pn, vol


def pressure_solver(Vtarget, phase, displacement, history, pressure, vol_tol, method='root_scalar'):
    """
    Solve for pressure using scipy optimization methods.
    
    Args:
        method: 'brentq', 'root_scalar', 'minimize_scalar', or 'secant' (original method)
    """
    
    def objective_function(pn):
        """Objective function: difference between target and computed volume"""
        pressure.assign(pn)
        displacement.solve()
        unew = displacement.get()
        history.update(unew)
        
        VK = compute_fracture_volume(phase.get_old(), unew)
        return Vtarget - VK
    
    def volume_function(pn):
        """Volume function for root finding"""
        pressure.assign(pn)
        displacement.solve()
        unew = displacement.get()
        history.update(unew)
        
        return compute_fracture_volume(phase.get_old(), unew)
    
    # Initial pressure estimate
    pn_initial = float(pressure)
    
    try:

        # General root finding with different methods
        result = root_scalar(
            lambda p: volume_function(p) - Vtarget,
            x0=pn_initial,
            x1=pn_initial * 1.1,  # Second guess for secant method
            method='secant',
            xtol=vol_tol,
            maxiter=50
        )
        pn = result.root
        iterations = result.iterations
            

        # Final check
        final_error = abs(objective_function(pn)) / (abs(Vtarget) if abs(Vtarget) > 1e-15 else 1.0)

        if final_error > vol_tol:
            logger.warning(f"Method {method} converged with error {final_error:.2e} > {vol_tol:.2e}")
            return -1, pn
        
        logger.debug(f"Pressure solver ({method}) converged in {iterations} iterations")
        return iterations, pn
        
    except Exception as e:
        logger.error(f"Pressure solver ({method}) failed: {e}")
        return -1, pn_initial