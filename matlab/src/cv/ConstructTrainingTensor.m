function Ttrain = ConstructTrainingTensor(T, trainIdx)

toNan = setdiff(find(~isnan(T)), trainIdx);
Ttrain = T;
Ttrain(toNan) = NaN;

end