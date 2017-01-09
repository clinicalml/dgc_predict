x = [0 1 0 1 0 1]; 
y = [1 0 0 1 1 0]; 

d2 = CosDistMulti(x,y,2);
assert(norm(d2 - [1 0 1]) < 1e-3);

d3 = CosDistMulti(x,y,3);
assert(norm(d3 - [1 0.5]) < 1e-3);