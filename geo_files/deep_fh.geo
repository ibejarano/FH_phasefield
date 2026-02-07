Mesh.Algorithm = 8; 
gridsize = 5e-2;
ref_gridsize = 2.5e-3;
L = 0.5;
H = 1;
dy = ref_gridsize*4;
//
Point(1) = { -L , H/2  ,0.0, gridsize};
Point(2) = { -L , -H/2  ,0.0, gridsize};
Point(3) = {  L , -H/2  ,0.0, gridsize};
Point(4) = {  L , H/2  ,0.0, gridsize};
// cuadrado donde estará la fractura
Point(5) = { -L , dy ,0.0, ref_gridsize};
Point(6) = { -L , -dy ,0.0, ref_gridsize};
Point(7) = { L , -dy ,0.0, ref_gridsize};
Point(8) = { L , dy ,0.0, ref_gridsize};
//
Line(1) = {1,5};
Line(2) = {6,2};
Line(3) = {2,3};
Line(4) = {3,7};
Line(5) = {4,1};
Line(10) = {5,8};
Line(11) = {8,7};
Line(12) = {7,6};
Line(13) = {6,5};
Line(14) = {8,4};
//
Curve Loop(1) = {1,10,14,5};
Curve Loop(2) = {-12,-11,-10,-13};
Curve Loop(3) = {2,3,4,12};
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
//
// Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};