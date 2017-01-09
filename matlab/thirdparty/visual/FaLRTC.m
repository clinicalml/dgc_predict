%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Low Rank Tensor Completion (FaLRTC)
% Time: 03/11/2012
% Reference: "Tensor Completion for Estimating Missing Values 
% in Visual Data", PAMI, 2012.
% min_{X} : \Psi(X) 
% s.t.         : X_\Omega = M_\Omega
% \Psi(X) = max_{Z_{i(i)} <= 1}: <X, \sum_i Y_i> - 0.5 \mu_i \|Y_i\|_F^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% code 
function [Y, errList, testErr] = FaLRTC(M, Omega, Mtest, Omega_test, alpha, ...
    mu, L, C, maxIter, epsilon, X)

% initialization
if nargin < 11
    X= M;
    X(logical(1-Omega)) = mean(M(Omega));
end

assert(all(isnan(M(Omega_test))));


Y = X;
Z = X;
B = 0;

N = ndims(M);
dim = size(M);
Gx = zeros(dim);
errList = zeros(maxIter, 1);
testErr = zeros(maxIter, 1);
normTestErr = NaN;
Lmax = 10*sum(1./mu);
tmp = zeros(1, N);
for i = 1:N
    tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.3;
end
P = 1.15;
flatNum = 15;
slope = (tmp - mu) / (1-(maxIter-flatNum)^(-P));
offset = (mu*(maxIter-flatNum)^P - tmp) / ((maxIter-flatNum)^P-1); 

mu0 = mu;
for k = 1:maxIter
    if k > 1
         fprintf('FaLRTC: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    
    % update mu
    mu = max(slope / (k^P) +offset, mu0);
    
    a2m = alpha.^2 ./ mu;
    ma = mu ./ alpha;
    
    Ylast = Y;
    %%  test L
    %L = L*C;
    while true
        b = (1+sqrt(1+4*L*B)) / (2*L);
        X = b/(B+b) * Z + B/(B+b) * Ylast;
        
        % compute f'(x) namely "Gx" and f(x) namely "fx"
        Gx = Gx * 0;
        fx = 0;
        for i = 1 : N
            [temp, sigma2] = Truncate(Unfold(X, dim, i), ma(i));
            temp = Fold(temp, dim, i);
            Gx = Gx + a2m(i) * temp;
            fx = fx + a2m(i)*(sum(sigma2) - sum(max(sqrt(sigma2)-ma(i), 0).^2));
        end


        % compute f(Ytest) namely fy
        %Y_model = X - Gx / L; % before setting Gx equal to 0 on the known elements of the tensor
        Gx(Omega) = 0;
        Y = X - Gx / L;
        fy = 0;
        for i = 1 : N
            [sigma] = SingularValue(Unfold(Y, dim, i));
            fy = fy + a2m(i)*(sum(sigma.^2) - sum(max(sigma-ma(i), 0).^2));
        end
        if (fx - fy)*L < sum(Gx(:).^2)
            if L > Lmax
                fprintf('FaLRTC: iterations = %d   difference=%f\n Exceeded the Maximum Lipschitiz Constant\n\n', k, errList(k-1));
                errList = errList(1:k);
                return;
            end
            L = L/C;
        else
             break;
        end
    end
    
    errList(k) =  norm(Y(:)-Ylast(:));
    testErr(k) = ComputeError(Mtest, Y(Omega_test));
    if abs(norm(Mtest)- 1) < 1e-3
        normTestErr = ComputeError(Mtest, Y(Omega_test)/norm(Y(Omega_test)));
    end
    fprintf('FaLRTC: iteration %d, test error = %0.4f, after normalization: %0.4f\n', k, testErr(k), normTestErr);
    if errList(k) < epsilon
        break;
    end
    
    %% update Z, Y, and B
    Z = Z - b*Gx;
    B = B+b;
    
    if ~isreal(Y)
        disp('some elements of Y are not real!')
    end
end

errList = errList(1:k);
fprintf('FaLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
