
InitRand(3);
T = randn(3,10,3);
T(2,:,1) = NaN;
T(1,:,2) = NaN;

clear args;
args.alpha = 1/2;
args.knnd.K = 2;
args.knnc.K = 2;
T_est1 = CompleteTensorKNNDrugCell(T, args);

T1 = CompleteTensorKNNDrug(T, args.knnd);
T2 = CompleteTensorKNNCell(T, args.knnc);
T_est2 = (T1 + T2)/2;

assert(norm(T_est1(1,:,2) - T_est2(1,:,2)) < 1e-5);
assert(norm(T_est1(2,:,1) - T_est2(2,:,1)) < 1e-5);

clear args;
args.knnd.K = 2;
args.knnc.K = 2;
T_est1 = CompleteTensorKNNDrugCell(T, args);
disp('Warning: need to test KNNdc further (in particular, the selection of alpha');