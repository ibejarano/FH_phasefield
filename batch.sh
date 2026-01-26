#!/bin/bash
NAME=output_sym_close
GEO_NAME=sym_deep
mkdir $NAME
cp deep_fracture/$GEO_NAME.geo $NAME
cp kgd_2d_refactored.py $NAME
gmsh -2 ./$NAME/$GEO_NAME.geo -format msh2 -o ./$NAME/mesh.msh
dolfin-convert ./$NAME/mesh.msh ./$NAME/mesh.xml
python kgd_2d_refactored.py $NAME --symmetric
