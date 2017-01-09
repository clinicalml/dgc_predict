x = [1 0];
y = [0 1];
assert(CosDist(x,y) == 1);
z = [-1 0];
assert(CosDist(x,z) == 2);