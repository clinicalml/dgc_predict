% TENSORMIX_ADM - Computes the reconstruction of a partly 
%                 observed tensor via "Mixture" approach
%
% Syntax
%  function [X,Z,fval,gval]=tensormix_adm(X, I, yy, lambda, eta, tol, verbose)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function [X,Z,fval,gval]=tensormix_adm(X, I, yy, lambda, eta, tol, verbose)

if ~exist('tol','var')
  tol=1e-3;
end

if ~exist('eta','var')
  eta=1;
end

if ~exist('verbose','var')
  verbose=0;
end

sz=size(X);
nd=ndims(X);
ind=sub2ind(sz, I{:});
ind0=setdiff(1:prod(sz),ind);
m=size(yy,1);

Z=cell(1,nd);
V=cell(1,nd);
for jj=1:nd
  szj = [sz(jj), prod(sz)/sz(jj)];
  Z{jj} = zeros(szj);
  V{jj} = zeros(szj);
end

nsv=10*ones(1,nd);
kk=1;

alpha=yy;
A=zeros(sz); A(ind)=alpha;
while 1
  for jj=1:nd
    Ztmp = Z{jj}+eta*flatten(A,jj);
    [Z{jj},S{jj},nsv(jj)]=softth(Ztmp,eta,nsv(jj));
    V{jj}=(Ztmp-Z{jj})/eta;

    viol(jj)=norm(flatten(A,jj)-V{jj},'fro');
  end
 
  Zsum = zeros(sz);
  Vsum = zeros(sz);
  for jj=1:nd
    Zsum=Zsum+flatten_adj(Z{jj},sz,jj);
    Vsum=Vsum+flatten_adj(V{jj},sz,jj);
  end
  
  alpha = (yy-Zsum(ind)+eta*Vsum(ind))/(lambda+eta*nd);
  A(ind)=alpha;

   
  % Compute the objective
  if lambda>0
    fval(kk)=0.5*sum((Zsum(ind)-yy).^2)/lambda;
  else
    fval(kk)=0;
  end
  
  gval(kk)=norm(lambda*alpha - yy + Zsum(ind));

  for jj=1:nd
    fval(kk)=fval(kk)+sum(S{jj});
  end
  
  % Compute the dual objective
% $$$   fact=1;
% $$$   for jj=1:nd
% $$$     fact=min(fact, 1/norm(flatten(A,jj)));
% $$$   end
% $$$   aa = alpha*fact;
% $$$   
% $$$   if kk>1
% $$$     dval(kk)=max(dval(kk-1),-0.5*lambda*sum(aa.^2)+aa'*yy);
% $$$   else
% $$$     dval(kk)=-inf;
% $$$   end
  
  if verbose 
    fprintf('[%d] fval=%g gval=%g viol=%s\n', kk, fval(kk), ...
          gval(kk), printvec(viol));
%    fprintf('[%d] fval=%g dval=%g fact=%g\n', kk, fval(kk), ...
%          dval(kk), fact);
  end
%  if kk>1 && 1-dval(kk)/fval(kk)<tol
  if kk>1 && max(viol)<tol && gval(kk)<tol
   break;
  end
  
  if kk>1000
    break;
  end

  kk=kk+1;
end

fprintf('[%d] fval=%g gval=%g viol=%s\n', kk, fval(kk), ...
        gval(kk), printvec(viol));

%fprintf('[%d] fval=%g dval=%g fact=%g\n', kk, fval(kk), ...
%        dval(kk), fact);

X=zeros(sz);
for jj=1:nd
  X = X + flatten_adj(Z{jj},sz,jj);
  Z{jj}=Z{jj}*nd;
end
