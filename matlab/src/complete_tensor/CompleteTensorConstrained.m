function [T_model, out] = CompleteTensorConstrained(T, args)

sz = size(T);
ind = find(~isnan(T));
[I,J,K] = ind2sub(sz, ind);
yy = T(ind);

[T_model,~,~,out.f,out.g] = tensorconst_adm(zeros(sz),...
                         {I,J,K}, yy, args.lambda, args.eta, args.tol);
                 

end
