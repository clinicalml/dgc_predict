InitRand(3);

T = randn(3,10,3);
T(2,:,1) = NaN;
T(1,:,2) = NaN;


%M = Unfold(T, size(T), 3)';
%C = corr(M, 'rows', 'pairwise');

C = nan(size(T,3),size(T,3),size(T,1)); 
for d = 1:size(T,1)
    C(:,:,d) = corr(squeeze(T(d,:,:))); 
end
C = mean(C, 3, 'omitnan');

args.K = 2;
T_est = CompleteTensorKNNCell(T, args);

w = C(2,[1 3]);
x = (w(1) * T(1,:,1) + w(2) * T(1,:,3)) / sum(w);
assert(norm(x - T_est(1,:,2)) < 1e-5);

w = C(1, [2 3]);
x = (w(1) * T(2,:,2) + w(2) * T(2,:,3)) / sum(w);
assert(norm(x - T_est(2,:,1)) < 1e-5);