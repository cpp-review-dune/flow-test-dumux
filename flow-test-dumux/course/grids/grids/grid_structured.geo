X = 60;
Y = 40;
Res = 1;

Point(1) = {0,0,0,Res};
Point(2) = {X,0,0,Res};
Point(3) = {X,Y,0,Res};
Point(4) = {0,Y,0,Res};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};

Plane Surface(10) = 9;

numCellsX = 24;
numCellsY = 16;

Transfinite Line{5} = numCellsX + 1;
Transfinite Line{6} = numCellsY + 1;
Transfinite Line{7} = numCellsX + 1;
Transfinite Line{8} = numCellsY + 1;

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

