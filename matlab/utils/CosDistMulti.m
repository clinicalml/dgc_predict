function [d,ny,nx] = CosDistMulti(x, y, len)
% This assumes that x and y are each a concatenation of vectors, each of
% length len, and hence computes the cosine distance for each set of
% vectors

% For my purposes, len would be nGenes

assert(mod(length(x),len)==0);
assert(length(x) == length(y));

n = length(x) / len;

for i = 1:n
    idx_start = (i-1)*len + 1;
    idx_end = i*len;
    [d(i), ny(i), nx(i)] = CosDist(x(idx_start:idx_end), y(idx_start:idx_end));
end


end
