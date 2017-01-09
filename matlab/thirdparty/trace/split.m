% SPLIT - Splits a string (or a cell array of strings) into a cell array
%
% Syntax
%  function B = split(A, sep)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function B = split(A, sep)

if ischar(A)
  A = {A};
end

B = cell(size(A));

for i=1:prod(size(A))
  s = findstr(A{i},sep);
  c = 1;
  B{i} = [];
  for j=1:length(s)
    B{i} = [B{i}, {A{i}(c:s(j)-1)}];
    c = s(j)+1;
  end
  if (c<length(A{i}))
    B{i} = [B{i}, {A{i}(c:end)}];
  end
end

if prod(size(A))==1
  B = B{1};
end
