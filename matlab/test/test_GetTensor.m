

[T,pertIds, geneIds, cellIds] = GetTensor('small_g50');

assert(size(T,1) == length(pertIds));
assert(size(T,2) == length(geneIds));
assert(size(T,3) == length(cellIds));
