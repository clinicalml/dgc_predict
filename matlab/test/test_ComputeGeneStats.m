

T = ones(3,3,2);
T(1,:,1) = NaN;

[m,v] = ComputeGeneStats(T);

assert(all(m == 1));
assert(all(v == 0));
