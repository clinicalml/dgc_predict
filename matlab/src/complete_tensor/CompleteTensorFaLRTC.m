function [T_model, out] = CompleteTensorFaLRTC(T, args)
    Omega = ~isnan(T);
    [T_model, ~, out.test_err] = FaLRTC(T, Omega, args.Mtest, args.Omega_test, args.alpha, ...
        args.mu, args.L0, args.C, args.maxIter,args.epsilon); 
end
