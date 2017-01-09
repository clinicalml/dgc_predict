T = GetTensor('d6_24hr');
clear args;
idx_obs = find(~isnan(T));

args = GetArgs('unfoldmc', [], [], [], []);
T1 = CompleteTensorUnfoldMC(T, args);

% check that original elements remain unchanged
assert(norm(T1(idx_obs) - T(idx_obs)) == 0);

% check that entire tensor is completed
assert(isempty(find(isnan(T1))));

%% check that heldout predictions are positively correlated with truth
% holdout 100 signatures
[T_remain, sigs, pairs, T_removed] = RemoveSigsRandom(T, 100); 

% complete tensor
T_complete = CompleteTensorUnfoldMC(T_remain, args);

% look at correlation
idx_removed = find(~isnan(T_removed));
assert(corr(T_complete(idx_removed), T_removed(idx_removed)) > 0.4);
