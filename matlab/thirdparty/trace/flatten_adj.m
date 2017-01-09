% FLATTEN_ADJ - mode-k folding of a matrix X into a tensor
%
% Syntax
%  function X=flatten_adj(X,sz,k)
%
% See also
%  FLATTEN
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function X=flatten_adj(X,sz,jj)
nd=length(sz);
sz=sz([jj:nd, 1:jj-1]);
X=reshape(X,sz);
if jj~=1
  X=permute(X,[nd-jj+2:nd 1:nd-jj+1]);
end
