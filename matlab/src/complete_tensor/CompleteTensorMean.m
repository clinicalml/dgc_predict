function [T_model, out] = CompleteTensorMean(T, args)

A = nanmean(T,3);
T_mean = repmat(A,[1 1 size(T,3)]);
T_model = T;
T_model(isnan(T)) = T_mean(isnan(T));
out = [];

end