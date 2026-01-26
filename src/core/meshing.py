import os
import shutil
import subprocess
import logging

logger = logging.getLogger(__name__)

def generate_mesh(geo_path: str, output_dir: str, mesh_name: str = "mesh") -> str:
    """
    Generates a FEniCS XML mesh from a Gmsh .geo file.
    
    1. Copies .geo file to output_dir
    2. Runs gmsh -2 to generate .msh
    3. Runs dolfin-convert to generate .xml
    
    Args:
        geo_path: Path to the source .geo file
        output_dir: Directory where the mesh will be generated
        mesh_name: Base name for the mesh files (default: "mesh")
        
    Returns:
        str: Absolute path to the generated XML file
    """
    if not os.path.exists(geo_path):
        raise FileNotFoundError(f"Geometry file not found: {geo_path}")
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Paths
    geo_dest = os.path.join(output_dir, f"{mesh_name}.geo")
    msh_path = os.path.join(output_dir, f"{mesh_name}.msh")
    xml_path = os.path.join(output_dir, f"{mesh_name}.xml")
    
    # 1. Copy geo file
    logger.info(f"Copying geometry {geo_path} -> {geo_dest}")
    shutil.copy(geo_path, geo_dest)
    
    # 2. Run Gmsh
    logger.info("Running Gmsh...")
    try:
        cmd_gmsh = ["gmsh", "-2", geo_dest, "-format", "msh2", "-o", msh_path]
        subprocess.run(cmd_gmsh, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logger.error(f"Gmsh failed: {e.stderr.decode()}")
        raise RuntimeError("Gmsh generation failed")
        
    # 3. Run dolfin-convert
    logger.info("Running dolfin-convert...")
    try:
        cmd_convert = ["dolfin-convert", msh_path, xml_path]
        subprocess.run(cmd_convert, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logger.error(f"dolfin-convert failed: {e.stderr.decode()}")
        raise RuntimeError("Mesh conversion failed")
        
    logger.info(f"Mesh generated successfully: {xml_path}")
    return xml_path
