import subprocess
from math import sin, cos
from dolfin import *
import numpy as np, math
from src.core.physics import elastic_energy_functional

def dem_KI_axisym(u, a, mu, nu, tip_z=0.0,
                  rho_min_frac=0.02, rho_max_frac=0.10, npts=10,
                  half_model=True, p=None):
    """
    DEM (displacement extrapolation) para K_I en penny-shaped axisimétrico.

    Parámetros
    ----------
    u : Function (Vector, 2D)   Desplazamiento resuelto (u_r, u_z)
    a : float                   Radio de la grieta
    mu: float                   Módulo cortante (G=μ)
    nu: float                   Poisson
    tip_z : float               Coordenada z del frente (típicamente 0.0)
    rho_min_frac, rho_max_frac  Rango de muestreo en [ρ_min, ρ_max] = [f_min*a, f_max*a]
    npts : int                  Nº de puntos en el ajuste lineal
    half_model : bool           Si usaste media placa (simetría), toma Δu=2*u_z
    p : float or None           Si lo das, devuelve también F = K_I/(p*sqrt(pi*a))

    Devuelve
    --------
    out : dict con keys {'a0','KI','F'(opcional),'samples'}
    """
    mesh = u.function_space().mesh()
    # Constantes
    kappa = 3.0 - 4.0*nu
    eps_z = 1e-8*max(1.0, float(a))  # pequeño desplazamiento por encima de z=0

    # Construir muestreo de distancias ρ y puntos r = a - ρ
    rhos = np.linspace(rho_min_frac*a, rho_max_frac*a, npts)
    rs   = a - rhos

    uz_vals = []
    for r_i in rs:
        try:
            uz = u(Point(float(r_i), float(tip_z + eps_z)))[1]
        except Exception:
            # Si cae justo fuera de un elemento, mover un pelo hacia adentro
            uz = u(Point(float(r_i - 1e-12), float(tip_z + eps_z)))[1]
        uz_vals.append(uz)

    uz_vals = np.array(uz_vals)
    # Apertura (Δu) en media placa ~ 2 * u_z (normal a la grieta)
    if half_model:
        delta_u = 2.0*uz_vals
    else:
        # en modelo completo necesitarías u_z^+ - u_z^- (no implementado aquí)
        delta_u = 2.0*uz_vals

    y = delta_u / np.sqrt(rhos)                # y = Δu / sqrt(ρ)
    # Ajuste lineal y(ρ) = a0 + b ρ
    coeffs = np.polyfit(rhos, y, 1)
    b, a0 = coeffs[0], coeffs[1]

    # K_I = (μ * sqrt(2π)/(κ+1)) * a0
    KI = (mu * math.sqrt(2.0*math.pi) / (kappa + 1.0)) * a0

    out = {"a0": float(a0), "KI": float(KI),
           "samples": {"rho": rhos, "delta_u": delta_u, "y": y, "fit": {"a0": a0, "b": b}}}

    if p is not None and a > 0.0:
        F = KI / (p * math.sqrt(math.pi * a))
        out["F"] = float(F)
    return out

def reemplazar_H(ruta_geo: str,  Lcrack: float, H_nuevo: float = None, beta_nuevo: float = None) -> None:
    # Leer todas las líneas
    with open(ruta_geo, 'r', encoding='utf-8') as f:
        lineas = f.readlines()

    # Escribir de nuevo, reemplazando solo la línea que empieza con "H ="
    with open(ruta_geo, 'w', encoding='utf-8') as f:
        for linea in lineas:
            if linea.lstrip().startswith("H ="):
                # Conservamos el mismo formato: "H = <valor>;"
                if H_nuevo is not None:
                    f.write(f"H = {H_nuevo};\n")

            elif linea.lstrip().startswith("beta_angle ="):
                if beta_nuevo is not None:
                    f.write(f"beta_angle = {beta_nuevo};\n")
            elif linea.lstrip().startswith("Lcrack ="):
                f.write(f"Lcrack = {Lcrack};\n")
            else:
                f.write(linea)

def run_gmsh(mesh_name: str, Lcrack=None , H_prof = None, beta = None , mesh=True):
    if mesh:
        reemplazar_H(f"{mesh_name}.geo", Lcrack, H_prof, beta)
        cmd_mallado = ["gmsh", "-2", f"{mesh_name}.geo", "-format", "msh2", "-o", f"{mesh_name}.msh"] 
        subprocess.run(cmd_mallado, check=True)
        cmd_xml_transformer = ["dolfin-convert" ,f"{mesh_name}.msh", f"{mesh_name}.xml"]
        subprocess.run(cmd_xml_transformer, check=True)
    
    mesh = Mesh(f"{mesh_name}.xml")
    boundaries = MeshFunction("size_t", mesh, f"{mesh_name}_facet_region.xml")
    return mesh, boundaries


def run_shallow_case(H_prof: float, 
                     p1: float, 
                     pxx: float,
                     Lcrack: float,
                     E: float,
                     nu:float,
                     beta: float = 0,
                     geo_name="shallow.geo", 
                     save_vtu=False,
                     tol= 1e-5,
                     mesh=True):
    from dolfin import ds
    np.seterr(divide='raise')

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack=Lcrack, H_prof=H_prof, beta=beta, mesh=mesh)

    check_H_sup = np.max(mesh.coordinates()[:, 1])

    print(f"Simulando caso: {check_H_sup:.4f} (malla) {H_prof:.4f} (dato) mts. profundidad")

    V = VectorFunctionSpace(mesh, "P", 2)

    bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
    bcs = [bc_bottom]


    mu = E / (2 * (1 + nu))

    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
    # lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = elastic_energy_functional(u, v, lmbda, mu)

    n = FacetNormal(mesh)
    lateral_compression = Constant((pxx, 0.0))
    L_form  = - dot(p1*n, v)*ds(10) - dot(p1*n, v)*ds(11)
    L_form += dot(lateral_compression, v)*ds(1) - dot(lateral_compression, v)*ds(3)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)

    # INICIO DEL POSTPROCESO
    kappa = 3 - 4 *nu # Plane - strain 
    # kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
    rel_KI = (p1 * np.sqrt(np.pi * Lcrack))

    npoints = 500

    beta_rad = np.deg2rad(beta)
    xs = np.linspace(Lcrack*0.75*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
    ys = np.linspace(Lcrack*0.75*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    KI_calc = np.zeros(npoints)
    KII_calc = np.zeros(npoints)

    r_crack = np.zeros(npoints)
    d_offset = tol
    if save_vtu:
        file = File('shallow_lateral.pvd')
        file << u_sol
    for i, (x, y) in enumerate(zip(xs, ys)):
        U = u_sol(x-d_offset, y+d_offset)
        Un = (U[1] * cos(beta_rad) - U[0] * sin(beta_rad))
        Ut = (U[0] * cos(beta_rad) + U[1] * sin(beta_rad))

        uplus_res[i] = [Ut, Un]

        V = u_sol(x+d_offset, y-d_offset)
        Vn = (V[1] * cos(beta_rad) - V[0] * sin(beta_rad))
        Vt = (V[0] * cos(beta_rad) + V[1] * sin(beta_rad))

        uminus_res[i] = [Vt, Vn]

        try:
            r_crack[i] = Lcrack -  np.sqrt(x**2 + y**2)
            KI_calc[i] = factor * abs(Un - Vn)/np.sqrt(r_crack[i])
            KII_calc[i] = factor * abs(Ut - Vt)/np.sqrt(r_crack[i])
        except FloatingPointError:
            r_crack[i] = 0
            KI_calc[i] = 0
            KII_calc[i] = 0

    calc_KI = np.max(KI_calc)/rel_KI
    calc_KII = np.max(KII_calc)/rel_KI
    #K_calcs[j] = [calc_KI, calc_KII]
    us = u_sol(0, check_H_sup)[1]
    up = u_sol(0, d_offset)[1]
    um = u_sol(0, -d_offset)[1]


    return [calc_KI, calc_KII], [us, up, um]

def run_shallow_case_symm(H_prof: float, 
                     p1: float, 
                     pxx: float,
                     Lcrack: float,
                     E: float,
                     nu:float,
                     beta: float = 0,
                     geo_name="shallow.geo", 
                     save_vtu=False,
                     tol= 1e-5,
                     plane_stess=False,
                     mesh=True):
    from dolfin import ds
    np.seterr(divide='raise')

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack=Lcrack, H_prof=H_prof, beta=beta, mesh=mesh)

    check_H_sup = np.max(mesh.coordinates()[:, 1])

    print(f"Simulando caso: {check_H_sup:.4f} (malla) {H_prof:.4f} (dato) mts. profundidad")

    V = VectorFunctionSpace(mesh, "P", 2)

    bc_axis_y = DirichletBC(V.sub(0), Constant(0), boundaries, 1) # 1 es el tag de la cara izquierda
    bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
    bcs = [bc_bottom, bc_axis_y]


    mu = E / (2 * (1 + nu))

    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
    if plane_stess:
        lmbda = 2 * mu * nu / (1 - nu)

    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    espesor = 0.02
    a = espesor * elastic_energy_functional(u, v, lmbda, mu)

    n = FacetNormal(mesh)
    lateral_compression = Constant((pxx, 0.0))
    L_form  = - dot(p1*n, v)*ds(10) - dot(p1*n, v)*ds(11)
    L_form -= dot(lateral_compression, v)*ds(3)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)

    # INICIO DEL POSTPROCESO
    kappa = 3 - 4 *nu # Plane - strain
    if plane_stess:
        kappa = (3 - nu)/(1+nu)
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
    rel_KI = (p1 * np.sqrt(np.pi * Lcrack))

    npoints = 800

    beta_rad = np.deg2rad(beta)
    xs = np.linspace(Lcrack*0.9*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
    ys = np.linspace(Lcrack*0.9*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    KI_calc = np.zeros(npoints)
    KII_calc = np.zeros(npoints)

    r_crack = np.zeros(npoints)
    d_offset = tol
    if save_vtu:
        file = File('shallow_lateral.pvd')
        file << u_sol
    for i, (x, y) in enumerate(zip(xs, ys)):
        U = u_sol(x-d_offset, y+d_offset)
        Un = (U[1] * cos(beta_rad) - U[0] * sin(beta_rad))
        Ut = (U[0] * cos(beta_rad) + U[1] * sin(beta_rad))

        uplus_res[i] = [Ut, Un]

        V = u_sol(x+d_offset, y-d_offset)
        Vn = (V[1] * cos(beta_rad) - V[0] * sin(beta_rad))
        Vt = (V[0] * cos(beta_rad) + V[1] * sin(beta_rad))

        uminus_res[i] = [Vt, Vn]

        try:
            r_crack[i] = Lcrack -  np.sqrt(x**2 + y**2)
            KI_calc[i] = factor * abs(Un - Vn)/np.sqrt(r_crack[i])
            KII_calc[i] = factor * abs(Ut - Vt)/np.sqrt(r_crack[i])
        except FloatingPointError:
            r_crack[i] = 0
            KI_calc[i] = 0
            KII_calc[i] = 0

    calc_KI = np.max(KI_calc)/rel_KI
    calc_KII = np.max(KII_calc)/rel_KI
    #K_calcs[j] = [calc_KI, calc_KII]
    us = u_sol(0, check_H_sup)[1]
    up = u_sol(0, d_offset)[1]
    um = u_sol(0, -d_offset)[1]


    return [calc_KI, calc_KII], [us, up, um]


def run_shallow_axisym(H_prof: float, 
                     p1: float, 
                     pxx: float,
                     Lcrack: float,
                     E: float,
                     nu:float,
                     beta: float = 0,
                     geo_name="shallow.geo", 
                     save_vtu=False,
                     tol= 1e-5,
                     mesh=True):
    from dolfin import ds
    import ufl

    pi = 3.14159
    np.seterr(divide='raise')

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack=Lcrack, H_prof=H_prof, mesh=mesh)

    check_H_sup = np.max(mesh.coordinates()[:, 1])

    print(f"Simulando caso: {check_H_sup:.4f} (malla) {H_prof:.4f} (dato) mts. profundidad")

    # TODO: Por el momento func de interp orden 1 para testear
    V = VectorFunctionSpace(mesh, "P", 2)
    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_sol = Function(V, name="u")

    # Material isotropico 
    E  = Constant(E)
    nu = Constant(nu)
    mu = E/(2.0*(1.0+nu))
    lmbda = E*nu/((1.0+nu)*(1.0-2.0*nu))
    # coordenada radial
    x = SpatialCoordinate(mesh)
    r = x[0]
    eps_r = Constant(1e-12)          # "regularización" suave
    r_safe = ufl.sqrt(r*r + eps_r*eps_r)   # continuo y evita r=0 en denominador
    I3 = Identity(3)

    # Deformación axisimétrica (3x3) construida desde (u_r,u_z)
    def eps_axisym(w):
        # w es vectorial: (u_r, u_z)
        G = ufl.grad(w)              # 2x2: [[∂u_r/∂r, ∂u_r/∂z],
                                    #       [∂u_z/∂r, ∂u_z/∂z]]
        e_rr = G[0, 0]
        e_zz = G[1, 1]
        e_rz = 0.5*(G[0, 1] + G[1, 0])
        e_tt = w[0] / r_safe         # ε_θθ = u_r / r

        E = ufl.as_tensor([[e_rr, 0.0,  e_rz],
                        [0.0,  e_tt, 0.0 ],
                        [e_rz, 0.0,  e_zz]])
        return E

    def sigma_axisym(w):
        E = eps_axisym(w)
        trE = ufl.tr(E)
        return 2.0*mu*E + lmbda*trE*I3

    # Cargas de volumen/contorno (opcionales)
    b = Constant((0.0, 0.0))  # cuerpo [N/m^3] (ya 3D)
    t = Constant((0.0, p1))  # tracción axisimétrica [Pa]

    # Peso axisimétrico: 2π r
    weight = 2.0*pi*r_safe

    a = ufl.inner(sigma_axisym(u), eps_axisym(v)) * weight * dx
    L = dot(t, v)* weight * ds(10) - ufl.dot(t, v) * weight * ds(11) +  (ufl.dot(b, v) * weight * dx)
    axis = CompiledSubDomain("near(x[0], 0.0)")

    bc_axis = DirichletBC(V.sub(0), Constant(0.0), axis)
    bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
    bcs = [bc_bottom, bc_axis]

    solve(a == L, u_sol, bcs, solver_parameters={"linear_solver": "mumps"})

    out = dem_KI_axisym(u_sol, Lcrack, mu, nu, half_model=False, p=p1)

    if save_vtu:
        file = File('axisym.pvd')
        file << u_sol

    return out


def run_left_notch(
    p1: float, 
    Lcrack: float,
    E: float,
    nu:float,
    geo_name="caso_2.geo", 
    save_vtu=False,
    tol= 1e-5,
    mesh=True   
):
    

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack, mesh=mesh)

    V = VectorFunctionSpace(mesh, "P", 1)


    def midpoint(x, on_boundary):
        return near(x[0], 0.0) and (x[0] > Lcrack) 

    bc_bottom = DirichletBC(V.sub(0), Constant(0), boundaries, 2)
    bc_top = DirichletBC(V.sub(0), Constant(0), boundaries, 4)

    bcs = [bc_bottom, bc_top]

    mu = E / (2 * (1 + nu))

    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
    #lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

    ds = Measure("ds", subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = elastic_energy_functional(u, v, lmbda, mu)


    upper_traction = Constant((0.0, p1))
    L_form  = dot(upper_traction, v)*ds(4) - dot(upper_traction, v)*ds(2)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)


    npoints = 200

    xs = np.linspace(Lcrack*0.6, Lcrack*0.999, npoints)
    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    d_offset = tol

    for i, x in enumerate(xs):
        uplus_res[i] = u_sol(x, d_offset)
        uminus_res[i] = u_sol(x, -d_offset)

    dU = np.abs(uplus_res[:, 1] - uminus_res[:, 1]) # No TOCAR
    dV = np.abs(uplus_res[:, 0] - uminus_res[:, 0])

    kappa = 3 - 4 *nu # Plane - strain 
    kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)

    r = Lcrack-xs

    KI_est = np.zeros(npoints)
    KII_est = np.zeros(npoints)

    for i, r_x in enumerate(r):
        KI_est[i] = factor * np.sqrt(1/r_x) * dU[i]
        KII_est[i] = factor * np.sqrt(1/r_x) * dV[i]

    if save_vtu:
        file = File('caso2.pvd')
        file << u_sol

    return r, KI_est, KII_est