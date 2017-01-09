function [T, sigs, idx_removed] = RemoveSigsChosen(T, pairs)
% DESCRIPTION: This function removes signatures specified by 'pairs' and
% outputs the resulting tensor.
%
% INPUTS:
%  T            Tensor 
%  pairs        n x 2 list of drug-cell pairs (indexes into T)
%
% OUTPUTS: 
%  T            Tensor with selected signatures removed (set to NaN) 
%  sigs         Matrix of removed signatures (n x numGenes)
%  idx_removed  Indices of removed signatures

assert(size(pairs,2) == 2);

numGenes  = size(T,2);
n = size(pairs,1);
sigs = NaN(n, numGenes);

idx_removed = [];

for p = 1:n
    i = pairs(p,1);
    j = pairs(p,2);
    sigs(p,:) = T(i,:,j);
    T(i,:,j) = NaN;
    idx_removed = [idx_removed sub2ind(size(T), repmat(i,1,numGenes), 1:numGenes, ...
        repmat(j, 1, numGenes))];
end

end
 
 