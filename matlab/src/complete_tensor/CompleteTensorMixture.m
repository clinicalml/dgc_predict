function [T_model, out] = CompleteTensorMixture(T, args)
    sz = size(T);
    ind = find(~isnan(T));
    [I,J,K] = ind2sub(sz, ind);
    yy = T(ind);
    warning('need to double check that T_est is not set equal the input values');
    [T_model,~,out.f,out.g] = tensormix_adm(zeros(sz), {I,J,K}, yy, ...
        args.lambda, args.eta, args.tol);
    
end
