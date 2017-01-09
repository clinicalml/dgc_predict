
T = GetTensor('small_g50');

[test_data_vec, test_idx, T_sub, d_vec] = SplitTensorDenser(T);

%%% check size of outputs

assert(length(d_vec) == length(T_sub));

ten_perc_num_sigs = round(size(T,1) * size(T,3) * .1);
assert(length(test_idx) == ten_perc_num_sigs * size(T,2));

for d = 1:length(T_sub)
   assert(isequal(size(T_sub{d}), size(T))); 
end

%%% check test data

assert(isequal(T(test_idx), test_data_vec));

B = zeros(size(T));
B(test_idx) = 1;

A1 = B(:,1,:);
A2 = B(:,10,:);
assert(isequal(A1,A2));
assert(length(find(A1)) == ten_perc_num_sigs);

%%% check density of T_sub elements
n = [];
for d = 1:length(T_sub)
    assert(ComputeDensity(T_sub{d})*100 - d_vec(d) < 1);
    n(d) = NumSigs(T_sub{d});
end
diff = n(2:end) - n(1:end-1);
assert(all(diff == ten_perc_num_sigs));