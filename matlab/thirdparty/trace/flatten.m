% FLATTEN - mode-k unfolding of a tensor
%
% Syntax
%  function Z=flatten(X,k)
% 
% See also
%  FLATTEN_ADJ
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function Z=flatten(X,jj)

nd=ndims(X);

if jj~=1
  X=permute(X,[jj:nd, 1:jj-1]);
end

Z=X(:,:);