function [sigs, I, J] = RandSigs(T, n, replace)

%InitRand();

if ~exist('replace')
    replace = false;
end

A = squeeze(~isnan(T(:,1,:)));
indA = find(A);

samples = randsample(indA, n, replace);

[I,J] = ind2sub(size(A), samples);
sigs = NaN(n, size(T,2));
for i = 1:n
   sigs(i,:) = T(I(i),:,J(i)); 
end


end