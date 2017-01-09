function idx = MapPairs2Tensor(pairs, size3)

assert(size(pairs,2) == 2);

size2 = size3([1,3]);

I2 = sub2ind(size2, pairs(:,1), pairs(:,2));

idx = MapMatrixInd2Tensor(I2, size2, size3);

end