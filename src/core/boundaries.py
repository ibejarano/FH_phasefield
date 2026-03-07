from dolfin import DirichletBC, SubDomain, Constant, near, MeshFunction, Mesh
from typing import List, Tuple, Optional

class CrackDomain(SubDomain):
    def __init__(self, center: List[float], l0: float, w0: float):
        super().__init__()
        self.center = center
        self.l0 = l0
        self.w0 = w0

    def inside(self, x, on_boundary):
        return abs(x[0] - self.center[0]) <= self.l0 and abs(x[1] - self.center[1]) <= self.w0

class NoCrackZone(SubDomain):
    """Marks all points OUTSIDE a horizontal band around the crack plane.
    Used to enforce phi=0 (no damage) outside the band."""
    def __init__(self, crack_y: float, band_half_width: float):
        super().__init__()
        self.crack_y = crack_y
        self.hw = band_half_width

    def inside(self, x, on_boundary):
        return abs(x[1] - self.crack_y) > self.hw

def setup_boundary_conditions(
    phase_field, 
    displacement_field, 
    l_init: float,
    h_elem: float,
    crack_center: List[float] = [0.0, 0.0],
    upper_face_free: bool = False,
    symmetric: bool = False,
    axisymmetric: bool = False,
    phi_band_width: float = 0.0
) -> Tuple[List[DirichletBC], List[DirichletBC]]:
    """
    Configura las condiciones de borde para el problema de fractura hidráulica.
    
    Args:
        phase_field: Instancia del campo de fase.
        displacement_field: Instancia del campo de desplazamiento.
        l_init: Longitud inicial de la grieta.
        h_elem: Tamaño característico del elemento (para el ancho de la grieta inicial).
        crack_center: Centro de la grieta inicial [x, y].
        upper_face_free: Si es True, la cara superior es libre (Tracción). Si es False, u=0.
        symmetric: Si es True, aplica simetría en x=0 (u_x=0).
        axisymmetric: Si es True, aplica condiciones para problema axisimétrico.
    
    Returns:
        Tuple[List[DirichletBC], List[DirichletBC]]: Listas de BCs para [desplazamiento, fase].
    """
    
    V_phi = phase_field.V
    V_u = displacement_field.V
    mesh = phase_field.mesh
    
    coords = mesh.coordinates()
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    
    bcs_u = []
    bcs_phi = []

    # --- 1. Condiciones Mecánicas (Desplazamiento) ---
    
    if axisymmetric:
        # ============ AXISYMMETRIC CASE ============
        # Coordinates: (r, z) where r = x[0], z = x[1]
        
        # a. Eje r=0: u_r = 0 (evita singularidad 1/r)
        def axis_r0(x, on_boundary):
            return on_boundary and near(x[0], x_min)
        bc_axis = DirichletBC(V_u.sub(0), Constant(0.0), axis_r0)
        bcs_u.append(bc_axis)
                
        # c. Bottom z=z_min: u = 0 (clamped)
        def bottom_side(x, on_boundary):
            return on_boundary and near(x[1], y_min)
        bc_bottom = DirichletBC(V_u.sub(1), Constant(0.0), bottom_side)
        bcs_u.append(bc_bottom)
        
        # d. Top z=z_max: far field (clamped if not free)
        if not upper_face_free:
            def upper_side(x, on_boundary):
                return on_boundary and near(x[1], y_max)
            bc_upper = DirichletBC(V_u, Constant((0.0, 0.0)), upper_side)
            bcs_u.append(bc_upper)
    
    elif symmetric:
        # ============ SYMMETRIC PLANE STRAIN ============
        # Symmetry plane at x=0
        
        # a. Symmetry x=0: u_x = 0
        def symmetry_plane(x, on_boundary):
            return on_boundary and near(x[0], x_min)
        bc_symmetry = DirichletBC(V_u.sub(0), Constant(0.0), symmetry_plane)
        bcs_u.append(bc_symmetry)
        
        # b. Bottom y=y_min: u = 0
        def bottom_side(x, on_boundary):
            return on_boundary and near(x[1], y_min)
        bc_bottom = DirichletBC(V_u, Constant((0.0, 0.0)), bottom_side)
        bcs_u.append(bc_bottom)
        
        # c. Top (if not free)
        if not upper_face_free:
            def upper_side(x, on_boundary):
                return on_boundary and near(x[1], y_max)
            bc_upper = DirichletBC(V_u, Constant((0.0, 0.0)), upper_side)
            bcs_u.append(bc_upper)
    
    else:
        # ============ FULL DOMAIN ============
        
        # a. Bottom: clamped
        def bottom_side(x, on_boundary):
            return on_boundary and near(x[1], y_min)
        bc_bottom = DirichletBC(V_u.sub(1), Constant(0.0), bottom_side)
        bcs_u.append(bc_bottom)
        
        # b. Top (if not free)
        if not upper_face_free:
            def upper_side(x, on_boundary):
                return on_boundary and near(x[1], y_max)
            bc_upper = DirichletBC(V_u, Constant((0.0, 0.0)), upper_side)
            bcs_u.append(bc_upper)
        
        # c. Left/Right: roller (u_x = 0)
        def left_side(x, on_boundary):
            return on_boundary and near(x[0], x_min)
        def right_side(x, on_boundary):
            return on_boundary and near(x[0], x_max)
        bc_left = DirichletBC(V_u.sub(0), Constant(0.0), left_side)
        bc_right = DirichletBC(V_u.sub(0), Constant(0.0), right_side)
        bcs_u.append(bc_left)
        bcs_u.append(bc_right)

    # --- 2. Condiciones de Campo de Fase (Daño) ---
    
    # Grieta Inicial (Dirichlet phi=1 en la zona definida)
    crack_subdomain = CrackDomain(crack_center, l_init, h_elem)
    bc_phi_crack = DirichletBC(V_phi, Constant(1.0), crack_subdomain)
    bcs_phi.append(bc_phi_crack)
    
    # Restricción geométrica: phi=0 fuera de una banda horizontal
    if phi_band_width > 0.0:
        no_crack = NoCrackZone(crack_center[1], phi_band_width)
        bc_phi_band = DirichletBC(V_phi, Constant(0.0), no_crack)
        bcs_phi.append(bc_phi_band)
    
    return bcs_u, bcs_phi

def create_boundary_markers(mesh: Mesh, left_id: int = 10, right_id: int = 20) -> MeshFunction:
    """
    Crea marcadores de frontera para identificar lados izquierdo/derecho.
    Útil para aplicar tracciones o cargas Neumann.
    """
    markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    markers.set_all(0)

    class LeftBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].min())
        
    class RightBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].max())

    LeftBoundary().mark(markers, left_id)
    RightBoundary().mark(markers, right_id)
    
    return markers
