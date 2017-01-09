function [T_model, out] = CompleteTensorUnfoldMC(T, args)

X = Unfold(T, size(T), args.unfold_dim);

idx = find(~isnan(X));
nElts = numel(X);
y = X(idx);
M = opRestriction(nElts, idx);
M_model = NonCVX_MC(y, M, size(X), args.p);
T_model = Fold(M_model, size(T), args.unfold_dim);
out = [];

end
