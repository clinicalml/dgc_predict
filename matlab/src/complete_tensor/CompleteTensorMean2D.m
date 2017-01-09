function [T_model, out] = CompleteTensorMean2D(T, args)

alpha = args.alpha;

A = nanmean(T,3);
B = squeeze(nanmean(T,1));
Ta = repmat(A,[1 1 size(T,3)]);
Tb = repmat(B,[1 1 size(T,1)]);
Tb = permute(Tb, [3 1 2]); 

T_mean = (1-alpha)*Ta + alpha * Tb;
T_model = T;
T_model(isnan(T)) = T_mean(isnan(T));
out = [];

end