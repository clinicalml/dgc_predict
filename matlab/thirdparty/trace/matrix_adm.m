% MATRIX_ADM - Computes reconstruction of a partly observed matrix
%
% Syntax
%  function [X,Z,Y,fval,gval]=matrix_adm(X, I, Bv, lambda, eta, tol, verbose)
%
% See also
%  TENSOR_AS_MATRIX
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function [X,Z,Y,fval,gval]=matrix_adm(X, I, Bv, lambda, eta, tol, verbose)

if ~exist('tol','var')
  tol=1e-3;
end

if ~exist('verbose','var')
  verbose=0;
end

sz=size(X);

m=length(Bv);


Z=zeros(size(X));
Y=zeros(size(X));

B=zeros(sz);
ind=sub2ind(sz,I{:});
B(ind)=Bv;

nsv=10;
kk=1;
while 1
  if lambda>0
    X = (B/lambda+eta*Z-Y)./((B~=0)/lambda+eta);
  else
    X=Z-Y/eta;
    X(ind)=Bv;
  end
  
  Z0=Z;
  [Z,ss,nsv]=softth(X+Y/eta,1/eta,nsv);

  Y=Y+eta*(X-Z);
  
  viol = norm(X(:)-Z(:));
  

  fval(kk)=sum(svd(X));
  if lambda>0
    fval(kk)=fval(kk)+0.5*sum((X(ind)-Bv).^2)/lambda;
  end

  gval(kk)=eta*norm(Z(:)-Z0(:)); %norm(G(:));

  if verbose
    fprintf('[%d] fval=%g gval=%g viol=%g\n', kk, fval(kk), gval(kk), ...
            viol);
  end
  
  
  if gval(kk)<tol && viol<tol
    break;
  end
  
  if kk>2000;
    break;
  end

  kk=kk+1;
end

fprintf('[%d] fval=%g gval=%g viol=%g\n', kk, fval(kk), gval(kk), ...
        viol);
