
InitRand(3);
T = randn(3,10,3);
T(2,:,1) = NaN;
T(1,:,2) = NaN;

clear args;
args.K = 'foo';
T_est = CompleteTensorMean(T, args);

w = [1 1];
x = (w(1) * T(1,:,1) + w(2) * T(1,:,3)) / sum(w);
assert(norm(x - T_est(1,:,2)) < 1e-5);

w = [1 1];
x = (w(1) * T(2,:,2) + w(2) * T(2,:,3)) / sum(w);
assert(norm(x - T_est(2,:,1)) < 1e-5);