% TENSORCONST_ADM - Computes the reconstruction of a partly
%                   observed tensor via "Constraint" approach
%
% Syntax
%  function [X,Z,Y,fval,gval] = tensorconst_adm(X, I, Bv, lambda, eta, tol, verbose)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function [X,Z,Y,fval,gval] = tensorconst_adm(X, I, Bv, lambda, eta, tol, verbose)

if ~exist('tol','var')
  tol=1e-3;
end

if ~exist('verbose','var')
  verbose=0;
end

sz=size(X);
nd=ndims(X);
m =length(I{1});

if nd~=length(I)
  error('Number of dimensions mismatch.');
end

if m~=length(Bv)
  error('Number of samples mismatch.');
end

Z=cell(1,nd);
Y=cell(1,nd);
S=cell(1,nd);

for jj=1:nd
  szj = [sz(jj), prod(sz)/sz(jj)];
  Y{jj} = zeros(szj);
  Z{jj} = zeros(szj);
end

B=zeros(sz);
ind=sub2ind(sz, I{:});
B(ind)=Bv;

nsv=10*ones(1,nd);
kk=1;
while 1
  X1 = zeros(size(X));
  for jj=1:nd
    X1 = X1 - flatten_adj(Y{jj}-eta*Z{jj},sz,jj);
  end
  
  if lambda>0
    X1(ind) = X1(ind) + Bv/lambda;
    X=X1./((B~=0)/lambda + nd*eta);
  else
    X=X1/(eta*nd);
    X(ind)=Bv;
  end
  

  % Check derivative
% $$$   D=zeros(size(X));
% $$$   D(ind)=X(ind)-Bv;
% $$$   for jj=1:nd
% $$$     D=D+eta*flatten_adj(Y{jj}/eta+flatten(X,jj)-Z{jj},sz,jj);
% $$$   end
% $$$   fprintf('gnorm=%g\n',norm(D(:)));
  
  for jj=1:nd
    [Z{jj},S{jj},nsv(jj)] = softth(Y{jj}/eta+flatten(X,jj),1/eta,nsv(jj));
    
    % Check derivative
    % fprintf('max[%d]=%g\n',jj,max(svd(eta*(Z{jj}-flatten(X,jj)-Y{jj}/eta))));
  end

  for jj=1:nd
    V=flatten(X,jj)-Z{jj};
    Y{jj}=Y{jj}+eta*V;
    viol(jj)=norm(V(:));
  end
  
  
  % Compute the objective
  G=zeros(size(X));
  fval(kk)=0;
  for jj=1:nd
    fval(kk)=fval(kk)+sum(svd(flatten(X,jj)));
    G = G + flatten_adj(Y{jj},sz,jj);
  end
  if lambda>0
    fval(kk)=fval(kk)+0.5*sum((X(ind)-Bv).^2)/lambda;
    G(ind)=G(ind)+(X(ind)-Bv)/lambda;
  else
    G(ind)=0;
  end
  
  gval(kk)=norm(G(:));

  if verbose
    fprintf('k=%d fval=%g gnorm=%g viol=%s\n',...
            kk, fval(kk), gval(kk), printvec(viol));
  end
  
  if kk>1 && max(viol)<tol && gval(kk)<tol
    break;
  end
  
  if kk>1000
    break;
  end
  
  
  kk=kk+1;
end

fprintf('k=%d fval=%g gnorm=%g viol=%s\n',...
        kk, fval(kk), gval(kk), printvec(viol));
  