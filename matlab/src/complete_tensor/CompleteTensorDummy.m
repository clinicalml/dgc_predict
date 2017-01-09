function [T_model, out] = CompleteTensorDummy(T, args)

T_model = T;
T_model(isnan(T)) = 1;
out = [];

end