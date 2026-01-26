Mesh.Algorithm = 8; 
gridsize = 5e-1;
ref_gridsize = 1e-3;
L = 10;
L_refinado = 0.2;
H = 10;
dy = ref_gridsize*5;
//
Point(1) = { 0 , H/2  ,0.0, gridsize};
Point(2) = { 0 , -H/2  ,0.0, gridsize};
Point(3) = {  L , -H/2  ,0.0, gridsize};
Point(4) = {  L , H/2  ,0.0, gridsize};
// cuadrado donde estará la fractura
Point(5) = { 0 , dy ,0.0, ref_gridsize};
Point(6) = { 0 , -dy ,0.0, ref_gridsize};
Point(7) = { L_refinado , -dy ,0.0, ref_gridsize};
Point(8) = { L_refinado , dy ,0.0, ref_gridsize};
Point(9) = { 0, 0 ,0.0, ref_gridsize};
Point(10) = { L_refinado , 0.0 ,0.0, ref_gridsize};

//
Line(1) = {1,5};
Line(2) = {6,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,1};
Line(10) = {5,8};
Line(11) = {8,10};
Line(12) = {7,6};
Line(13) = {6,9};
Line(14) = {9,10};
Line(15) = {9,5};
Line(16) = {10,7};
//
Curve Loop(1) = {1, 10, 11, 16, 12, 2, 3, 4, 5};
Curve Loop(2) = {-15, 14, -11, -10};
Curve Loop(3) = {-13, -12, -16, -14 };
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
//
// Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};
