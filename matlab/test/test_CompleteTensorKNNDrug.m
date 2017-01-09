InitRand(3);
T = randn(3,10,3);
T(2,:,1) = NaN;
T(1,:,2) = NaN;

C = nan(size(T,1),size(T,1),size(T,3)); 
for c = 1:size(T,3)
    C(:,:,c) = corr(squeeze(T(:,:,c)')); 
end
C = mean(C, 3, 'omitnan');

%M = Unfold(T, size(T), 1)';
%C = corr(M, 'rows', 'pairwise');

clear args;
args.K = 2;
T_est = CompleteTensorKNNDrug(T, args);

w = C(1,[2 3]);
x = (w(1) * T(2,:,2) + w(2) * T(3,:,2)) / sum(w);
assert(norm(x - T_est(1,:,2)) < 1e-5);

w = C(2, [1 3]);
x = (w(1) * T(1,:,1) + w(2) * T(3,:,1)) / sum(w);
assert(norm(x - T_est(2,:,1)) < 1e-5);

%T_est2 = CompleteTensorKNNDrugSave(T, args);
%assert(norm(T_est(:) - T_est2(:)) < 1e-10);