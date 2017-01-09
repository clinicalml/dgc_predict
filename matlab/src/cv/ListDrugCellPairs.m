function [I_split, J_split] = ListDrugCellPairs(T, K, writeToFile, tensor_name)
% DESCRIPTION: This function is used for the purpose of creating input
% files to TensorCV2 or TensorCV3, for the purposes of parallelization.  It
% lists out drug-cell index tuples, splitting all signatures into K files.

A = squeeze(~isnan(T(:,1,:)));
indA = find(A);

samples = randsample(indA, length(indA), false);

[I,J] = ind2sub(size(A), samples);

% add padding if necessary
n = length(I);
if mod(n,K) ~= 0
    r = floor(n/K);
    nPad = K*(r+1) - n;
    I = [I; nan(nPad, 1)];
    J = [J; nan(nPad, 1)];
end

I_split = reshape(I, K, length(I) / K);
J_split = reshape(J, K, length(J) / K);

if writeToFile
   base = ['../data/results/' tensor_name '/cv/drug_cell_pairs/'];
   if ~exist(base)
       mkdir(base)
   end
       
   for i = 1:K
      pairs = [I_split(i,:); J_split(i,:)]';
      file = sprintf('%sdrug_cell_pairs_%d_of_%d.txt', base, i, K); 
      dlmwrite(file, pairs, ' ');
   end
end

end
