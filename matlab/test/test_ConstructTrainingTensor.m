

T = GetTensor('small_g50');
nSplits = 2;
[testIdx, trainIdx] = SplitTensor(T,nSplits);

Ttrain = ConstructTrainingTensor(T, trainIdx{1});

assert(ComputeDensity(Ttrain) / ComputeDensity(T) - 0.9 < 1e-3);

idx = find(~isnan(Ttrain));
assert(all(Ttrain(idx) == T(idx)));



