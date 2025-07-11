//Works as intended!
// Top left square
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

// Top right square
Point(5) = {2, 0, 0, 1.0};
Point(6) = {2, 1, 0, 1.0};

// Bottom right square
Point(7) = {1, -1, 0, 1.0};
Point(8) = {2, -1, 0, 1.0};

// Lines for top left square
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Lines for top right square
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

// Lines for bottom right square
Line(8) = {2, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};

// Curve loops and surfaces
Curve Loop(1) = {1, 2, 3, 4};           // Top left
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, -2};           // Top right
Plane Surface(2) = {2};

Curve Loop(3) = {8, 9, 10, -5};          // Bottom right
Plane Surface(3) = {3};

// Structured mesh
Transfinite Line {1,2,3,4,5,6,7,8,9,10} = 11 Using Progression 1;
Transfinite Surface {1,2,3};
Recombine Surface {1,2,3}; // Optional: quadrangles

// Physical groups (optional)
Physical Surface("TopLeft") = {1};
Physical Surface("TopRight") = {2};
Physical Surface("BottomRight") = {3};
