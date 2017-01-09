function [T_model, out] = CompleteTensorSiLRTC(T, args)
    Omega = ~isnan(T);
    [T_model, ~, out.test_err] = SiLRTC(T, Omega, args.Ttest, ...
        args.Omega_test, args.alpha, args.beta, args.maxIter, args.epsilon); 
end
