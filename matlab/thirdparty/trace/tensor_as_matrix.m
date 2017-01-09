% TENSOR_AS_MATRIX - Computes the reconstruction of partly observed
%                    tensor via "As A Matrix" approach
%
% Syntax
%  function [X,Z,fval,gval]=tensor_as_matrix(X, I, Bv, eta, tol);
%
% See also
%  MATRIX_ADM
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function [X,Z,fval,gval]=tensor_as_matrix(X, I, Bv, eta, tol)

if ~exist('tol','var')
  tol=1e-3;
end

sz=size(X);
nd=ndims(X);

Z=cell(1,nd);
for ii=1:nd
  szp=[sz(ii:end) sz(1:ii-1)];
  Ip=[I(ii:end) I(1:ii-1)];
  J =sub2ind(szp(2:end), Ip{2:end});
  [Z1,Z{ii},Y,fval1,gval1]=matrix_adm(zeros(szp(1),prod(szp(2:end))),{Ip{1}, ...
                      J}, Bv, 0, eta, tol);
  fval(ii)=fval1(end);
  gval(ii)=gval1(end);
end

X=[];