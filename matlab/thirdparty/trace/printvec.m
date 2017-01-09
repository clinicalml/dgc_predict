% PRINTVEC - Prints a vector into a string
%
% Syntax
%  function str=printvec(vv,maxLength)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function str=printvec(vv,ml)

if ~exist('ml','var')
  ml=10;
end


vv=vv(:);

str = '[';

for ii=1:min(length(vv),ml)-1
  str = [str, sprintf('%g ', vv(ii))];
end

if length(vv)>ml
  str = [str, '...'];
end

str = [str, sprintf('%g]',vv(end))];
