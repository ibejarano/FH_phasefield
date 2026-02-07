// Simple small mesh for testing axisymmetric formulation
// Domain: r in [0, 1], z in [0, 1]
// multiplico todo por 5

Mesh.Algorithm = 5;
Mesh.RandomFactor = 1e-3;
Mesh.Smoothing = 3;

gridsize = 1;
ref_gridsize = 4e-3;

L = 20.0;  // radial extent
H = 20.0;  // axial extent
H_top = 1;

// Outer corners
Point(1) = { 0 , H_top  , 0.0, gridsize};
Point(2) = { 0 , -H/2 , 0.0, gridsize};
Point(3) = { L , -H/2 , 0.0, gridsize};
Point(4) = { L , H_top , 0.0, gridsize};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// definimos la curva

Point(5) = {0.005, 0.0, 0};
Point(6) = {0.4, 0.0, 0};
Point(7) = {0.76, 0.03, 0};
Point(8) = {1.23, 0.11, 0};
Point(9) = {1.95, 0.41, 0};
Point(10) = {2.7, 0.82, 0};

// Construccion de curvas de refinamiento
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
// Curve {5, 6} In Surface {1};

// 3. Configurar el refinamiento mediante Fields
// Campo de Distancia: mide la distancia a la línea 5
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8, 9};
Field[1].Sampling = 100;

// Campo Threshold: define el tamaño en función de la distancia
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = ref_gridsize;  // Tamaño sobre la curva
Field[2].SizeMax = gridsize;           // Tamaño lejos de la curva
Field[2].DistMin = 0.05;         // Distancia donde empieza a ensanchar
Field[2].DistMax = 3;         // Distancia donde ya alcanza el tamaño máximo

// 4. Establecer como campo de fondo
Background Field = 2;

// Opcional: Evitar que los elementos sean muy alargados
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
