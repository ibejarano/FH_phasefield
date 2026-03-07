// Simple small mesh for testing axisymmetric formulation
// Domain: r in [0, 1], z in [0, 1]
// multiplico todo por 5

Mesh.Algorithm = 5;
Mesh.RandomFactor = 1e-3;
Mesh.Smoothing = 3;

scale_case = 10;

gridsize = 1;
ref_gridsize = 3e-3/scale_case;

L = 20.0;  // radial extent
H = 20.0;  // axial extent
H_top = 1;
L_frac = 3/scale_case;

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

Point(5) = {0.0001/scale_case, 0.0, 0};
Point(6) = {L_frac, 0.0, 0};

// Construccion de curvas de refinamiento
Line(5) = {5, 6};
Curve {5} In Surface {1};

// 3. Configurar el refinamiento mediante Fields
// Campo de Distancia: mide la distancia a la línea 5
Field[1] = Distance;
Field[1].CurvesList = {5};
Field[1].Sampling = 200;

// Campo Threshold: define el tamaño en función de la distancia
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = ref_gridsize;  // Tamaño sobre la curva
Field[2].SizeMax = gridsize;           // Tamaño lejos de la curva
Field[2].DistMin = 0.02/scale_case;         // Distancia donde empieza a ensanchar
Field[2].DistMax = 3;         // Distancia donde ya alcanza el tamaño máximo

// 4. Establecer como campo de fondo
Background Field = 2;

// Opcional: Evitar que los elementos sean muy alargados
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
