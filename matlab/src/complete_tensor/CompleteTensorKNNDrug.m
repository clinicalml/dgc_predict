function [T_model, out] = CompleteTensorKNNDrug(T, args)
    
T = permute(T, [3 2 1]);
[T_model, out] = CompleteTensorKNNCell(T, args);
T_model = permute(T_model, [3 2 1]);

end
