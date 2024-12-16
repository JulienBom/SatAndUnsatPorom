lc = 0.06;

tot_length = 4;
tot_height = 4;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {tot_length,0.0,0.0,lc};
Point(3) = {tot_length,tot_height,0.0,lc};
Point(4) = {0.0,tot_height,0.0,lc};

Line(10) = {4,3};
Line(20) = {3,2};
Line(30) = {2,1};
Line(40) = {1,4};
Line Loop(50) = {10,20,30,40};
Plane Surface(60) = {50};
