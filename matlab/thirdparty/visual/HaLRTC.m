%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADM algorithm: tensor completion
% paper: Tensor completion for estimating missing values in visual data
% date: 05-22-2011
% min_X: \sum_i \alpha_i \|X_{i(i)}\|_*
% s.t.:  X_\Omega = T_\Omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X_est, errList, testErr] = HaLRTC(T, Omega, Ttest, Omega_test, alpha, ...
    beta, maxIter, epsilon, X)

assert(all(isnan(T(Omega_test))));

if nargin < 9
    X = T;
    X(logical(1-Omega)) = mean(T(Omega));
end

if length(Ttest) > 0 && length(Omega_test) > 0
    compute_test_error = true;
else
    compute_test_error = false;
end

errList = zeros(maxIter, 1);
testErr = zeros(maxIter, 1);
dim = size(T);
Y = cell(ndims(T), 1);
M = Y;

normT = norm(T(:));
for i = 1:ndims(T)
    Y{i} = X;
    M{i} = zeros(dim);
end

Msum = zeros(dim);
Ysum = zeros(dim);
for k = 1: maxIter
    if mod(k, 2) == 0
        if compute_test_error
        fprintf('HaLRTC: iterations = %d   difference=%f, test_err=%f, cos_err=%f, mean_norm_true=%f, mean_norm-est=%f\n',...
            k, errList(k-1), testErr(k-1), cosErr, meanNormTrue, meanNormEst);
        else
            fprintf('HaLRTC: iterations = %d   difference=%f, \n', k, errList(k-1));
        end
    end
    beta = beta * 1.05;
    
    % update Y
    Msum = 0*Msum;
    Ysum = 0*Ysum;
    for i = 1:ndims(T)
        Y{i} = Fold(Pro2TraceNorm(Unfold(X-M{i}/beta, dim, i), alpha(i)/beta), dim, i);
        Msum = Msum + M{i};
        Ysum = Ysum + Y{i};
    end
    
    % update X
    %X(logical(1-Omega)) = ( Msum(logical(1-Omega)) + beta*Ysum(logical(1-Omega)) ) / (ndims(T)*beta);
    lastX = X;
    X = (Msum + beta*Ysum) / (ndims(T)*beta);
    X_est = X;
    X(Omega) = T(Omega);
    
    % update M
    for i = 1:ndims(T)
        M{i} = M{i} + beta*(Y{i} - X);
    end
    
    % compute errors
    errList(k) = norm(X(:)-lastX(:)); 
    
    if compute_test_error
        testErr(k) = ComputeError(Ttest, X(Omega_test));
        [cc, n_est, n_true] = CosDistMulti(Ttest, X(Omega_test), size(T, 2));
        cosErr = mean(cc);
        meanNormEst = mean(n_est);
        meanNormTrue = mean(n_true);
    end
    
    if errList(k) < epsilon
        break;
    end
end

errList = errList(1:k);
fprintf('HaLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));

