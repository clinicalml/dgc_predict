function [T_model, out] = CompleteTensorConstrained(T, args)

sz = size(T);
ind = find(~isnan(T));
[I,J,K] = ind2sub(sz, ind);
yy = T(ind);

warning('need to check whether T_model is set to equal yy at the known indices');

[T_model,~,~,out.f,out.g] = tensorconst_adm(zeros(sz),...
                         {I,J,K}, yy, args.lambda, args.eta, args.tol);
                 

end