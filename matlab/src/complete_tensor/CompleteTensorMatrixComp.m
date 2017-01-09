function [T_model, out] = CompleteTensorMatrixComp(T, args)


X = squeeze(T(:,1,:));
IDX = find(~isnan(X));
nElts = numel(X);
T_model = T;

for i=1:size(T,2)
   s = sprintf('i=%d',i);
   disp(s);
   X = squeeze(T(:,i,:)); 
   y = X(IDX);
   M = opRestriction(nElts, IDX);
   T_model(:,i,:) = NonCVX_MC(y, M,size(X), args.p);
end

out = []; % I should output the cross-validated hyperparameter?

end