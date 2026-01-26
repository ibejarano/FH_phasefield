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
    symmetric: bool = False
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
    
    # a. Borde Inferior: Empotrado o deslizante según el caso (aquí fijo u=0)
    #    Para simulaciones KGD 'deep', solemos fijar u_y=0 o u=0 en bordes lejanos.
    #    Aquí se mantiene la lógica original: fijo completo en el fondo.
    def bottom_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].min())

    bc_bottom = DirichletBC(V_u, Constant((0.0, 0.0)), bottom_side)
    bcs_u.append(bc_bottom)

    # b. Borde Superior
    if not upper_face_free:
        def upper_side(x, on_boundary):
            return near(x[1], mesh.coordinates()[:, 1].max())
        
        bc_upper = DirichletBC(V_u, Constant((0.0, 0.0)), upper_side)
        bcs_u.append(bc_upper)

    # c. Simetría
    if symmetric:
        def symmetry_plane(x, on_boundary):
            return on_boundary and near(x[0], 0.0)
        
        # Restringe u_x = 0 en el plano de simetría (subespacio 0 del vector)
        bc_symmetry = DirichletBC(V_u.sub(0), Constant(0.0), symmetry_plane)
        bcs_u.append(bc_symmetry)

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
