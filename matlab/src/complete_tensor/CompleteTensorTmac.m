function [T_model, out] = CompleteTensorTmac(T, args)

data = T(~isnan(T));
known = find(~isnan(T));
Nway = size(T);

[X,Y,out] = TMac(data, known, Nway, args.estCoreNway, args.opts);   

T_model = zeros(Nway);
for n = 1:length(Nway)
    if out.alpha(n) > 0
        Tn = Fold(X{n}*Y{n},Nway,n);
        T_model = T_model + out.alpha(n) * Tn;
    end
end


end