% Demo

% Problem Size: n1, n2 - matrix dimensions, r - rank
n1 = 50; n2 = 60; r = 5;

% Generating matrix
X = randn(n1,r)*randn(r,n2);

% Vectrozizing the matrix for projection
x = X(:);

sizeX = size(X);

% Generating random sampling points
T = randperm(prod(sizeX));
IDX = T(1:round(0.4*prod(sizeX))); % 40% sampling

% Creating operator for selecting entries at the chosen random locations 
% Requires Sparco Toolbox 
% http://www.cs.ubc.ca/labs/scl/sparco/
M = opRestriction(prod(sizeX), IDX);

% Sampled data
y = M(x,1);

% For proposed algorithms

XRec = NonCVX_MC(y,M,sizeX,0.1); % Non Convex Algorithm for p = 0.1
norm(X-XRec,'fro')/norm(X,'fro') % Normalized Mean Squared Error

XRec = IHT_MC(y,M,sizeX); % Regularized Iterated Hard Thresholding (worst results)
norm(X-XRec,'fro')/norm(X,'fro') % Normalized Mean Squared Error

XRec = IST_MC(y,M,sizeX); % Regularized Iterated Soft Thresholding
norm(X-XRec,'fro')/norm(X,'fro') % Normalized Mean Squared Error

% For Comparison with other algorithms
% Need to download SVT (Singular Value Thresholding) Toolbox 
% http://svt.caltech.edu/

% The format for the SVT inputs is slightly different
% We keep the same variable names as used in SVT
Omega = IDX;
data = X(Omega);

% Parameters as suggested in SVT Toolbox
tau = 5*sqrt(n1*n2);
p  = length(IDX)/prod(sizeX);
delta = 1.2/p;
mu_final = .01; 
tol = 1e-3;

% Solving via Singular Value Thresholding algorithm
[U,S,V,numiter] = SVT([n1 n2],Omega,data,tau,delta);
XRec = U*S*V';
norm(X-XRec,'fro')/norm(X,'fro') % Normalized Mean Squared Error

% Solving via Fixed Point Continuation algorithm
[U,S,V,numiter] = FPC([sizeX(1) sizeX(2)],Omega,data,mu_final);
XRec = U*S*V';
norm(X-XRec,'fro')/norm(X,'fro') % Normalized Mean Squared Error