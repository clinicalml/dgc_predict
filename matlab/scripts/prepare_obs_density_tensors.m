[X,pertIds_full,geneIds,cellIds] = GetTensor('tsize/small/small');

nDrugs = size(X,1);

n = round(linspace(50,nDrugs,5));

for i = 1:(length(n)-1)
    for rep = 1:5
    idx = sort(randsample(nDrugs, n(i), false)); 
    T = X(idx,:,:);
    pertIds = pertIds_full(idx);
    filename = sprintf('%s/expr/tensor/tsize/small/obs_density/T%d_rep%d.mat', DataDir(), i, rep);
    save(filename, 'T', 'pertIds', 'geneIds', 'cellIds');
    end
end
