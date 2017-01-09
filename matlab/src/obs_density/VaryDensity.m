function ind_keep = VaryDensity(T, d)

ind_keep = {};
A = squeeze(~isnan(T(:,1,:)));
ind2 = find(A);
maxsigs = size(T,1)*size(T,3);

for i = 1:length(d)
    nsigs_to_keep = round(maxsigs*d(i));
    p = randperm(length(ind2), nsigs_to_keep);
    ind2_keep = ind2(p);
    ind_keep{i} = MapMatrixInd2Tensor(ind2_keep, size(A), size(T));
end

end
