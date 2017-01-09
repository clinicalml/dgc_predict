function [T_model, out] = CompleteTensorAsMatrix(T, args)

    sz = size(T);
    ind = find(~isnan(T));
    [I,J,K] = ind2sub(sz, ind);
    yy = T(ind);
    
    [~,Z,out.f,out.g] = tensor_as_matrix(zeros(sz), {I,J,K}, yy, args.eta, args.tol);
    
    for j = 1:3
        out.X{j} = flatten_adj(Z{j},sz,j);
    end
    
    args.alpha = args.alpha / sum(args.alpha);
    T_model = args.alpha(1)*out.X{1} + args.alpha(2)*out.X{2} + args.alpha(3)*out.X{3};
end
