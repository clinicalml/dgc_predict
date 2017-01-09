function X = IHT_MC(y,M,sizeX,err,x_initial,normfac,insweep,tol,decfac)

% Matrix Completion via Iterated Hard Thresholding
% min rank(X) subject to ||y - M(X)||_2<err

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries
% err - norm of the mismatch (default 0)
% x_initial - intial estimate of the vectorized form of X (defalut 0)
% normfac - eigenvalue of (M'M) (default should be 1 for masking operator)
% insweep - maximum number of internal sweeps for solving ||y - M(X)||_2 + lambda rank(X) (default 200)
% tol - tolerance (default 1e-4)
% decfac - decrease factor for cooling lambda

% Copyright (c) Angshul Majumdar 2010

if nargin < 4
    err = 1e-6;
end
if nargin < 5
    x_initial = zeros(prod(sizeX),1);
end
if nargin < 6
    normfac = 1;
end
if nargin < 7
    insweep = 200;
end
if nargin < 8
    tol = 1e-4;    
end
if nargin < 9
    decfac = 0.9;
end

alpha = 1.1*normfac;
x = x_initial;
lambdaInit = decfac*max(abs(M(y,2))); lambda = lambdaInit;
f_current = norm(y-M(x,1)) + lambda*norm(x,1);

while lambda > lambdaInit*tol
    for ins = 1:insweep
        f_previous = f_current;
        x = x + (1/alpha)*M(y - M(x,1),2);
        [U,S,V] = svd(reshape(x,sizeX),'econ');
        s = HardTh(diag(S),lambda/(2*alpha));
        S = diag(s);
        X = U*S*V';
        x = X(:);
        f_current = norm(y-M(x,1)) + lambda*norm(x,1);
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
    end
    if norm(y-M(x,1))<err
        break;
    end
    lambda = decfac*lambda;
end

    function  z = HardTh(s,thld)
        z = (abs(s) >= thld) .* s;
    end
end
