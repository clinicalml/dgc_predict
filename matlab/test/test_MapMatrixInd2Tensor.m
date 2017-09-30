

rng(123,'twister');
T = randn(4,3,2);
toNan = randperm(8, 5);

toNan3 = MapMatrixInd2Tensor(toNan, [4 2], [4 3 2]);
T(toNan3) = NaN;

for i = 1:3
   A = T(:,i,:);
   I = sort(find(isnan(A)));
   assert(all(I' == sort(toNan)));
end
