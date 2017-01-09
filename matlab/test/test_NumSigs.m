
T = GetTensor('T_test');
n = NumSigs(T);
assert(n == 36);

d = ComputeDensity(T);
nPossible = size(T,1) * size(T,3);
assert( abs(n/nPossible - d) < 1e-3);

perDrug = NumSigs(T, 'drug');
assert(all(perDrug' == [8 7 8 7 6]));

perCell = NumSigs(T, 'cell');
assert(all(perCell == [5 4 5 2 5 5 5 5]));