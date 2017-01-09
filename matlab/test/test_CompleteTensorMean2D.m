
InitRand(3);
T = randn(3,10,3);
T(2,:,1) = NaN;
T(1,:,2) = NaN;

clear args;
args.alpha = 0.5;
T_est = CompleteTensorMean2D(T, args);

x = (T(1,:,1) + T(1,:,3)+ T(2,:,2) + T(3,:,2)) / 4;
assert(norm(x - T_est(1,:,2)) < 1e-5);

x = (T(1,:,1) + T(3,:,1)+ T(2,:,2) + T(2,:,3)) / 4;
assert(norm(x - T_est(2,:,1)) < 1e-5);