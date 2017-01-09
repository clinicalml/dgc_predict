function [errRate] = ComputeErrRate(T1, T2)

x1 = T1(:);
x2 = T2(:);

idx1 = find(isnan(x1));
idx2 = find(isnan(x2));

assert(isequal(idx1, idx2));

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

errRate = norm(x1-x2) / norm(x1);

end