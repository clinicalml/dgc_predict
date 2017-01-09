function I3 = MapMatrixInd2Tensor(I2, size2, size3)
% DESCRIPTION: This function takes a vector of linear indices of some
% length n, corresponding to a matrix of dimensions 'size2' (specifically,
% the drug-cell matrix), and maps them to a vector of linear indices (of
% length n * numGenes) for all the corresponding elements of the tensor,
% extended along the gene dimension.
% 
% INPUTS:
% I2        vector of linear indices into drug-cell matrix
% size2     vector of length 2 specifying the dimensions of the drug-cell matrix
% size3     vector of length 3 specifying the dimensions of the tensor
% 
% OUTPUTS:
% I3        vector of linaer indices into the tensor

% map into a tensor that extends the second dimension
assert(size2(1) == size3(1));
assert(size2(2) == size3(3));

A = NaN * ones(size2);
A(I2) = 1;

B = repmat(A,[1 1 size3(2)]);
B = permute(B, [1 3 2]);
assert(all(size(B) == size3));

I3_v1 = find(~isnan(B));

[I,J,K] = ind2sub(size3, I3_v1);
M = [I J K];
M2 = sortrows(M,[1 3]);

I3 = sub2ind(size3, M2(:,1), M2(:,2), M2(:,3));

end