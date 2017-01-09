function [T_model, out] = CompleteTensorHaLRTC(T, args)
    Omega = ~isnan(T);
    [T_model, ~, out.test_err] = HaLRTC(T, Omega, args.Ttest, ...
        args.Omega_test, args.alpha, args.beta, args.maxIter, args.epsilon); 
end
