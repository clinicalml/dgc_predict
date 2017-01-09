T = reshape(1:24, 2, 3, 4);
dim = size(T);

M1 = Unfold(T, dim, 1);
M2 = Unfold(T, dim, 2);
M3 = Unfold(T, dim, 3);

T1 = Fold(M1, dim, 1);
T2 = Fold(M2, dim, 2);
T3 = Fold(M3, dim, 3);

assert(isequaln(T, T1));
assert(isequaln(T, T2));
assert(isequaln(T, T3));

T4 = Fold(M3, dim, 1);
assert(~isequaln(T, T4));