from dolfin import Mesh, MPI
import os

def check_mesh(path):
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return
    
    mesh = Mesh(path)
    hmin = mesh.hmin()
    hmax = mesh.hmax()
    num_cells = mesh.num_cells()
    
    # Parameters logic from Model
    relacion_hl = 4.0
    l_c = hmin * relacion_hl
    l_init = l_c * 2.5
    
    print(f"--- {path} ---")
    print(f"Num Cells: {num_cells}")
    print(f"hmin: {hmin:.6e}")
    print(f"hmax: {hmax:.6e}")
    print(f"l_c (sim): {l_c:.6e}")
    print(f"l_init (sim): {l_init:.6e}")
    print("")

check_mesh("output_full/mesh.xml")
check_mesh("output_sym/mesh.xml")
