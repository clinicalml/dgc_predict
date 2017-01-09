
rng('default');
rng(123);

T = GetTensor('small_g50');
nSplits = 2;
[testIdx, trainIdx] = SplitTensor(T,nSplits);
n = length(find(~isnan(T)));

for i = 1:nSplits
    
    t = sort(T([testIdx{i};trainIdx{i}]));

    if(i > 1)
        assert(all(oldt==t));
    end
    assert(all(~isnan(t)));
    assert(length(t) == n);
    
    oldt = t;

end

p = length(t) / length(find(~isnan(T)));
assert(1-p < 1e-3);

T_add = GetTensor('small_g50');
T2 = cat(1, T, T_add);

rng('default');
rng(123);

pertIdx = 1:size(T_add,1);

[testIdx2, trainIdx2] = SplitTensor(T2, nSplits, nSplits, pertIdx);
n2 = length(find(~isnan(T2)));

% the testIdx should have the same drug,cell coordinates
for i = 1:nSplits
    [I, J, K ] = ind2sub(size(T), testIdx{i});
    [I2,J2,K2] = ind2sub(size(T2), testIdx2{i});
    assert(isequal(I,I2));
    assert(isequal(J,J2));
    assert(isequal(K,K2));
    
    assert(length(trainIdx2{i}) == 3*length(trainIdx{i}));
    assert(length([trainIdx2{i}; testIdx2{i}]) == n2);
    assert(isequal(sort(union(trainIdx2{i},testIdx2{i})), find(~isnan(T2))));
end

% M = ConstructTrainingTensor(T2, testIdx2{1});
% B = ~isnan(M(1:3,1:6,1:9));
% binaryTensorVoxel(B);

