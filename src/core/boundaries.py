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

def setup_boundary_conditions(
    phase_field, 
    displacement_field, 
    l_init: float,
    h_elem: float,
    crack_center: List[float] = [0.0, 0.0],
    upper_face_free: bool = False,
    symmetric: bool = False,
    axisymmetric: bool = False
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
    
    Returns:
        Tuple[List[DirichletBC], List[DirichletBC]]: Listas de BCs para [desplazamiento, fase].
    """
    
    V_phi = phase_field.V
    V_u = displacement_field.V
    mesh = phase_field.mesh
    
    bcs_u = []
    bcs_phi = []

    # --- 1. Condiciones Mecánicas (Desplazamiento) ---
    
    # a. Borde Inferior (Bottom)
    # Case A: Full Domain (symmetric=False). Bottom is far-field. Fix u=0.
    # Case B: Symmetric Domain (symmetric=True). Bottom is Crack Plane (z=0).
    #         - Condition: u_y = 0 on ligament (x > l_init).
    #         - Condition: Free on crack (x < l_init).
    #         - Note: In Phase Field, clamping u_y=0 right at the tip might affect diffusion.
    #         - We will apply u_y=0 for x > l_init.
    
    if not symmetric:
        # Far field bottom -> Clamped or simple support
        def bottom_side(x, on_boundary):
            return near(x[1], mesh.coordinates()[:, 1].min())
        bc_bottom = DirichletBC(V_u, Constant((0.0, 0.0)), bottom_side)
        bcs_u.append(bc_bottom)
    else:
        # Symmetric Half -> Bottom is symmetry plane
        # 1. Ligament (x > l_init): u_y = 0
        def ligament_side(x, on_boundary):
            is_bottom = near(x[1], mesh.coordinates()[:, 1].min())
            is_ligament = x[0] > l_init # Strictly greater to allow tip opening? 
            # Or use >= but maybe l_init node needs to open if damage is there.
            # Using > l_init ensures the node exactly at l_init is NOT constrained.
            return on_boundary and is_bottom and is_ligament
        
        # Constrain u_y = 0 on ligament
        bc_ligament_y = DirichletBC(V_u.sub(1), Constant(0.0), ligament_side)
        bcs_u.append(bc_ligament_y)
        
        # 2. Prevent rigid body motion in X if needed? 
        # In Axisymmetric, x=0 is fixed (u_x=0).
        # In Plane Strain Symmetric, x=0 is usually symmetry plane too (u_x=0).
        # So X is constrained by the vertical boundaries.
        pass

    # b. Borde Superior
    if not upper_face_free:
        def upper_side(x, on_boundary):
            return near(x[1], mesh.coordinates()[:, 1].max())
        
        bc_upper = DirichletBC(V_u, Constant((0.0, 0.0)), upper_side)
        bcs_u.append(bc_upper)

    # c. Simetría Vertical / Axisimetría (x=0)
    if symmetric or axisymmetric:
        def axis_boundary(x, on_boundary):
            return on_boundary and near(x[0], 0.0)
        
        # Restringe u_x = 0 (Roller)
        # In Axi: u_r = 0. In Sym: u_x = 0. Consistent.
        bc_axis = DirichletBC(V_u.sub(0), Constant(0.0), axis_boundary)
        bcs_u.append(bc_axis)

    if not symmetric:
        def left_side(x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].min())
        def right_side(x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].max())
        bc_left = DirichletBC(V_u.sub(0), Constant(0.0), left_side)
        bc_right = DirichletBC(V_u.sub(0), Constant(0.0), right_side)
        bcs_u.append(bc_left)
        bcs_u.append(bc_right)

    # --- 2. Condiciones de Campo de Fase (Daño) ---
    
    # Grieta Inicial (Dirichlet phi=1 en la zona definida)
    crack_subdomain = CrackDomain(crack_center, l_init, h_elem)
    bc_phi_crack = DirichletBC(V_phi, Constant(1.0), crack_subdomain)
    bcs_phi.append(bc_phi_crack)
    
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
