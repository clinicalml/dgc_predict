
InitRand(123);
T1 = RandTensor([10 10 10], 3, 0);
T2 = T1 + randn(10, 10, 10)*0.2;
err = ComputeErrRate(T1, T2);

% err should be comparable to the 0.2 noise scaling, right?
assert(err > 0.1); 
assert(err < 0.3);
