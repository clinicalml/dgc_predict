function Ttrain = ConstructTrainingTensorFrom2D(T, trainIdx_2d)

A = squeeze(T(:,1,:));
trainIdx_3d = MapMatrixInd2Tensor(trainIdx_2d, size(A), size(T));
toNan = setdiff(find(~isnan(T)), trainIdx_3d);
Ttrain = T;
Ttrain(toNan) = NaN;

end