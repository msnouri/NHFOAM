// Circular basin, R = 1 m, depth 0.5 m: unstructured triangles extruded
// into a sigma mesh of prisms (vertical columns)
SetFactory("OpenCASCADE");

R = 1.0;
h = 0.5;
lc = 0.06;
nLayers = 8;

Disk(1) = {0, 0, -h, R, R};
Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;

out[] = Extrude {0, 0, h} { Surface{1}; Layers{nLayers}; Recombine; };

Physical Surface("surface") = {out[0]};
Physical Surface("bottom") = {1};
Physical Surface("wall") = {out[2]};
Physical Volume("water") = {out[1]};
