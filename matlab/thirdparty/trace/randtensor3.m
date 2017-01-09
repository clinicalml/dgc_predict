% RANDTENSOR3 - Randomly generates a 3-way low-rank tensor
%
% Syntax
%  function X=randtensor3(sz, dims)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function X=randtensor3(sz, dims)

nd=length(sz);

C=randn(dims);

U=cell(1,nd);
for jj=1:nd
  [U{jj},R]=qr(randn(sz(jj),dims(jj)),0);
end

X=kolda3(C,U{:});
