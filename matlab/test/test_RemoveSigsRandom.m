% setup
T = GetTensor('small_g50', false);
N = 1000;
InitRand();

% run function for both 'keep1' = true and false
[A_remaining, ~, ~, A_removed] = RemoveSigsRandom(T, N, true);
[B_remaining, ~, ~, B_removed] = RemoveSigsRandom(T, N, false);

% combine remaining and removed tensors and check that I get T back
A_remaining(isnan(A_remaining)) = 0;
A_removed(isnan(A_removed)) = 0;
A_full = A_remaining + A_removed;
T_test = T;
T_test(isnan(T)) = 0;
assert(isequaln(T_test, A_full));

B_remaining(isnan(B_remaining)) = 0;
B_removed(isnan(B_removed)) = 0;
B_full = B_remaining + B_removed;
T_test = T;
T_test(isnan(T)) = 0;
assert(isequaln(T_test, A_full));


% just rerun since I messed up A_remaining and A_removed
T = GetTensor('small_g50', false);
N = 1000;
InitRand();
[A_remaining, sigsA, pairsA, A_removed] = RemoveSigsRandom(T, N, true);
[B_remaining, sigsB, pairsB, B_removed] = RemoveSigsRandom(T, N, false);

% check number of (nonNan) elements of output tensors
a = length(find(~isnan(A_remaining(:))));
b = length(find(~isnan(B_remaining(:))));
assert(a == b);
assert(a == (NumSigs(T) - N) * size(T,2));

% check num sigs of output tensors
assert(NumSigs(A_remaining) == NumSigs(T) - N);
assert(NumSigs(B_remaining) == NumSigs(T) - N);

% check size of sig matrices
assert(isequal(size(sigsA), [N, size(T,2)]));
assert(isequal(size(sigsB), [N, size(T,2)]));

% check that I have no NaNs in sigs matrices
assert(sum(isnan(sigsA(:))) == 0);
assert(sum(isnan(sigsB(:))) == 0);

% check that all pairs are unique
assert(size(unique(pairsA, 'rows'), 1) == N);
assert(size(unique(pairsB, 'rows'), 1) == N);

% check that linear indexing gives same numbers
idxA = MapPairs2Tensor(pairsA, size(A_remaining));
assert(isequal(sort(A_removed(idxA)),sort(sigsA(:))));
