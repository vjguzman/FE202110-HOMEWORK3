// Gmsh project created on Wed May  5 23:32:26 2021
SetFactory("OpenCASCADE");

L    = 16.0;
H    = 4.0;
r    = 1.0;
Lext = 2.0;

Point(1) = {-L/2, -H/2, 0, mesh};
Point(2) = {L/2,  -H/2, 0, mesh};
Point(3) = {L/2,   H/2, 0, mesh};
Point(4) = {-L/2,  H/2, 0, mesh};


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};


//+
Circle(5) = {0, 0, 0, r, 0, 2*Pi};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
Extrude {-Lext, 0, 0} {
  Curve{4}; 
}
//+
Extrude {Lext, 0, 0} {
  Curve{2}; 
}

BooleanDifference{ Surface{1}; Delete; }{Surface{2}; Delete;}


//+
Physical Curve("Empotrado") = {8};
//+
Physical Curve("NaturalBoundary") = {11};
//+
Physical Surface("Placa") = {1};
//+
Physical Surface("Extremos") = {3, 4};
