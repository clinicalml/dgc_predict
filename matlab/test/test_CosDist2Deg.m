x = [0,1];
y = [1,0];
z = [-1,1];

d1 = CosDist(x,y);
d2 = CosDist(x,z);
d3 = CosDist(y,z);

tol = 1e-10;
assert(abs(CosDist2Deg(d1)-90) < tol);
assert(abs(CosDist2Deg(d2)-45) < tol);
assert(abs(CosDist2Deg(d3)-135) < tol);
