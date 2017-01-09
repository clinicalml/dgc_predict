function [U,V,Lambda,Theta,omega,pi, result] = phege(U0,V0,Lambda0,omega0,pi0,Theta0,P,G,R,lambdaU,lambdaV,deltaOmega,deltaPi,max_iter,tol)
U = U0;
V = V0;
Lambda = Lambda0;
omega = omega0;
pi = pi0;
Theta = Theta0;

l = eval_loss(U,V,Lambda,omega,pi,P,G,Theta,lambdaU,lambdaV,deltaOmega,deltaPi);

fprintf('iteration 0, loss = %f\n',l)

result(1,1) = l;

for i=1:max_iter
    W = U*Lambda*V';
    Theta = solve_Theta(W,R);
    %fprintf('iteration %d,updating Theta ...\n',i);
    
    Sigma = U*U';
    omega = solve_omega(P,Sigma,deltaOmega);
    %fprintf('iteration %d,updating omega ...\n',i);
    
    Sigma = V*V';
    pi = solve_omega(G,Sigma,deltaPi);
    %fprintf('iteration %d,updating pi ...\n',i);
    
    Lambda_Para = Setting_Lambda_Para(U,V,Theta,max_iter,tol);
    Lambda = solve_Lambda(Lambda,Lambda_Para);
    %fprintf('iteration %d,updating Lambda ...\n',i);
    
    U_Para = Setting_U_Para(V,Theta,Lambda,lambdaU,P,omega,max_iter,tol);
    U = solve_U(U,U_Para);
    %fprintf('iteration %d,updating U ...\n',i);
    
    V_Para = Setting_V_Para(U,Theta,Lambda,lambdaV,G,pi,max_iter,tol);
    V = solve_V(V,V_Para);
    %fprintf('iteration %d,updating V ...\n',i);
    
    l = eval_loss(U,V,Lambda,omega,pi,P,G,Theta,lambdaU,lambdaV,deltaOmega,deltaPi);
    fprintf('iteration %d,loss = %f ...\n',i,l);
    
    result(i+1, 1) = l;
    
end

%------evaluate total loss--------------------------------
function l = eval_loss(U,V,Lambda,omega,pi,P,G,Theta,lambdaU,lambdaV,deltaOmega,deltaPi)

A = Theta - U*Lambda*V';
l1 = sum(sum(A.*A));

numP = length(P);
a = zeros(numP,1);
for i=1:numP
    A = P{i}-U*U';
    a(i) = sum(sum(A.*A));
end
l2 = omega'*a;
r2 = sum(sum(omega.*omega));

numG = length(G);
a = zeros(numG,1);
for i=1:numG
    A = G{i}-V*V';
    a(i) = sum(sum(A.*A));
end
l3 = pi'*a;
r3 = sum(sum(pi.*pi));

l = l1+lambdaU*(l2+deltaOmega*r2)+lambdaV*(l3+deltaPi*r3);

%----------solving Theta-------------------------------------
function Theta = solve_Theta(W,R)
mask_nonzero = double(R>0);
mask_zero = double(R==0);
Theta = R.*mask_nonzero + W.*mask_zero;

%----------solving Lambda------------------------------------
function Lambda = solve_Lambda(Lambda0,Lambda_para)
[Lambda,grad,iter] = PGD(@eval_grad_Lambda,@eval_obj_Lambda,Lambda_para,Lambda0);

function obj_val = eval_obj_Lambda(Lambda,Lambda_para)
F = Lambda_para.Theta-Lambda_para.U*Lambda*Lambda_para.V';
obj_val = sum(sum(F.*F));

function grad = eval_grad_Lambda(Lambda,Lambda_para)
grad = 2*(Lambda_para.U'*Lambda_para.U*Lambda*(Lambda_para.V'*Lambda_para.V)-Lambda_para.U'*Lambda_para.Theta*Lambda_para.V);

function Lambda_Para = Setting_Lambda_Para(U,V,Theta,max_iter,tol)
Lambda_Para.U = U;
Lambda_Para.V = V;
Lambda_Para.Theta = Theta;
Lambda_Para.max_iter = max_iter;
Lambda_Para.tol = tol;

%-----------solving omega--------------------------------------
function omega = solve_omega(P,Sigma,delta)
% P is a cell array collecting all phenotype similarities
% solving pi is the same
num_similarity = length(P);
a = zeros(num_similarity,1);
for i = 1:num_similarity
    A = P{i}-Sigma;
    a(i) = sum(sum(A.*A));
end
omega = projsplx(a/(2*delta));

%-----------solving U-----------------------------------------
function U = solve_U(U0,U_para)
[U,grad,iter] = PGD(@eval_grad_U,@eval_obj_U,U_para,U0);

function obj_val = eval_obj_U(U,U_para)
F = U_para.Theta-U*U_para.Lambda*U_para.V';
f1 = sum(sum(F.*F));
num_similarity = length(U_para.P);
a = zeros(num_similarity,1);
UUt = U*U';
for i=1:num_similarity
    D = U_para.P{i}-UUt;
    a(i) = sum(sum(D.*D));
end
f2 = U_para.omega'*a;
obj_val = f1+U_para.lambda*f2;

function grad = eval_grad_U(U,U_para)
num_similarity = length(U_para.P);
Pt = U_para.P{1}*U_para.omega(1);
if num_similarity>1
    for i=2:num_similarity
        Pt = Pt+U_para.omega(i)*U_para.P{i};
    end
end
grad = 2*(U*U_para.Lambda*(U_para.V'*U_para.V)*U_para.Lambda'-U_para.Theta*U_para.V*U_para.Lambda'-U_para.lambda*Pt*U+2*U_para.lambda*(U*U')*U);

function U_Para = Setting_U_Para(V,Theta,Lambda,lambda,P,omega,max_iter,tol)
U_Para.V = V;
U_Para.Theta = Theta;
U_Para.Lambda = Lambda;
U_Para.lambda = lambda;
U_Para.P = P;
U_Para.max_iter = max_iter;
U_Para.tol = tol;
U_Para.omega = omega;

%-----------solving V-----------------------------------------
function V = solve_V(V0,V_para)
[V,grad,iter] = PGD(@eval_grad_V,@eval_obj_V,V_para,V0);

function obj_val = eval_obj_V(V,V_para)
F = V_para.Theta-V_para.U*V_para.Lambda*V';
f1 = sum(sum(F.*F));
num_similarity = length(V_para.G);
a = zeros(num_similarity,1);
VVt = V*V';
for i=1:num_similarity
    D = V_para.G{i}-VVt;
    a(i) = sum(sum(D.*D));
end
f2 = V_para.pi'*a;
obj_val = f1+V_para.lambda*f2;

function grad = eval_grad_V(V,V_para)
num_similarity = length(V_para.G);
Gt = V_para.G{1}*V_para.pi(1);
if num_similarity>1
    for i=2:num_similarity
        Gt = Gt+V_para.pi(i)*V_para.G{i};
    end
end
grad = 2*(V*V_para.Lambda'*(V_para.U'*V_para.U)*V_para.Lambda-V_para.Theta'*V_para.U*V_para.Lambda-V_para.lambda*Gt*V+2*V_para.lambda*(V*V')*V);

function V_Para = Setting_V_Para(U,Theta,Lambda,lambda,G,pi,max_iter,tol)
V_Para.U = U;
V_Para.Theta = Theta;
V_Para.Lambda = Lambda;
V_Para.lambda = lambda;
V_Para.G = G;
V_Para.max_iter = max_iter;
V_Para.tol = tol;
V_Para.pi = pi;

%----------general projected gradient-------------------------
function [X,grad,iter] = PGD(gradf,objeval,X_para,X0)

% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations

X = X0; 
maxiter = X_para.max_iter;
tol = X_para.tol;

alpha = 1; beta = 0.1;
obj_old = objeval(X,X_para);
for iter=1:maxiter,  
  grad = gradf(X,X_para);
  projgrad = norm(grad(grad < 0 | X >0));
  if projgrad < tol,
    break
  end

  % search step size 
  for inner_iter=1:20,
    Xn = max(X - alpha*grad, 0); d = Xn-X;
    gradd=sum(sum(grad.*d)); 
    obj_diff = objeval(Xn,X_para)-objeval(X,X_para);
    suff_decr = obj_diff - 0.01*gradd < 0;
    if inner_iter==1,
      decr_alpha = ~suff_decr; Xp = X;
    end
    if decr_alpha, 
      if suff_decr,
	X = Xn; break;
      else
	alpha = alpha * beta;
      end
    else
      if ~suff_decr | Xp == Xn,
	X = Xp; break;
      else
	alpha = alpha/beta; Xp = Xn;
      end
    end
  end
  obj_new = objeval(X,X_para);
  if obj_new == obj_old
      break;
  else
      obj_old = obj_new;
  end
  
end



%----------------general Euclidean Projection---------------------------
function w = ProjectOntoSimplex(v, b)
% PROJECTONTOSIMPLEX Projects point onto simplex of specified radius.
%
% w = ProjectOntoSimplex(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  sum(w) <= b, w >= 0.
%
%   That is, performs Euclidean projection of v to the positive simplex of
%   radius b.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)

if (b < 0)
  error('Radius of simplex is negative: %2.3f\n', b);
end
v = (v > 0) .* v;
u = sort(v,'descend');
sv = cumsum(u);
rho = find(u > (sv - b) ./ (1:length(u))', 1, 'last');
theta = max(0, (sv(rho) - b) / rho);
w = max(v - theta, 0);

