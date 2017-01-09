function [test_data_vec, test_idx, T_sub, d_vec] = ...
    SplitTensorDenser(T, d_min, d_incr, d_test)
% d stands for density (i.e. observation density)

%%% set defaults
if ~exist('d_min')
    d_min = 10;
end

if ~exist('d_incr')
    d_incr = 10;
end

if ~exist('d_test')
    d_test = 10;
end


%%% check inputs
assert(1 <= d_min && d_min < 100);
assert(1 <= d_incr && d_incr < 100);
assert(1 <= d_test && d_test < 100);


%%% remove extra sigs so that obs density is d_min plus a multiple of d_incr
d_tot = ComputeDensity(T)*100;
d_extra = mod(d_tot - d_min, d_incr);
n_extra = d2n(d_extra, size(T));
T = RemoveSigsRandom(T, n_extra, true);


%%% designate randomly selected test data
n_test = d2n(d_test, size(T));
[T_maxd, ~, test_pairs] = RemoveSigsRandom(T, n_test, true);
test_idx = MapPairs2Tensor(test_pairs, size(T));
test_data_vec = T(test_idx);


%%% incrementally remove more and more signatures randomly
d_max = round(ComputeDensity(T_maxd)*100);
d_vec = d_min:d_incr:d_max;
T_sub = {};
idx = length(d_vec);
T_sub{idx} = T_maxd;

n_incr = d2n(d_incr, size(T));
while (idx > 1)
    idx = idx - 1;
    T_sub{idx} = RemoveSigsRandom(T_sub{idx+1}, n_incr, true);
    
    % check that actual observation density is close to expected
    d_test = ComputeDensity(T_sub{idx})*100;
    assert(abs(d_test - d_vec(idx)) < 1); 
    
    % check that number of signatures per drug is at least 1
    assert(all(NumSigs(T_sub{idx}, 'drug') > 0))
end


end

function n = d2n(d, siz)
    n = round(siz(1)*siz(3)*d/100);
end


