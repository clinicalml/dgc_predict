% SOFTTH - Computes the proximity operator with respect to the
%          trace norm
%
% Syntax
%  function [vv,ss,nsv]=softth(vv, lambda, nsv, verbose);
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function [vv,ss,nsv]=softth(vv, lambda, nsv, verbose);

if ~exist('verbose','var')
  verbose=0;
end


sz=size(vv);
nsv=min(min(sz),nsv+1);

if verbose
  fprintf('sz=[%d %d]\n',sz(1), sz(2));
  fprintf('nsv=');
end
  
while 1
  if verbose
    fprintf('%d/',nsv);
  end
 [U,S,V]=pca(vv,min(min(sz),nsv),10);
 ss=diag(S);
  if min(ss)<lambda || nsv==min(sz)
    if verbose
      fprintf('min(ss)=%g\n',min(ss));
    end
    break;
  else
    nsv=min(min(sz),round(nsv*2));
  end
end

ix=find(ss>=lambda);
ss=ss(ix)-lambda;
vv = U(:,ix)*diag(ss)*V(:,ix)';

nsv=length(ix);