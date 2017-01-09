
InitRand(123);
T1 = randn(10,20,30);
pairs = [1 3; 2 5; 9 7];
[T2, sigs, idx_removed] = RemoveSigsChosen(T1, pairs);

assert(length(find(isnan(T2))) == 60);
assert(isequal(sigs, [T1(1,:,3); T1(2,:,5); T1(9,:,7)]));
assert(length(idx_removed) == 60);
assert(all(isnan(T2(idx_removed))));
idx_kept = setdiff(1:numel(T1), idx_removed);
assert(isempty(find(isnan(T1(idx_kept)), 1)));