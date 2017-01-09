% setup
T = GetTensor('small_g50', false);
N = 1000;
InitRand();

[T_remaining, sigs, pairs, T_removed] = RemoveSigsRandom(T, N, false);

sigs_test = zeros(N, size(T,2));
for p = 1:N
   sigs_test(p,:) = T(pairs(p,1), :, pairs(p,2));
end

assert(isequal(sigs, sigs_test));

idx = MapPairs2Tensor(pairs, size(T));
idx_test = find(~isnan(T_removed));

assert(isequal(sort(idx), sort(idx_test)));
assert(length(idx) == N*size(T,2));

A_full = squeeze(~isnan(T(:,1,:)));
A_remaining = squeeze(~isnan(T_remaining(:,1,:)));
A_removed = squeeze(~isnan(T_removed(:,1,:)));
idx_full = sort(find(A_full));
idx_remaining = find(A_remaining);
idx_removed = find(A_removed);
idx_union = sort(union(idx_remaining, idx_removed));

assert(isequal(idx_union, idx_full));
