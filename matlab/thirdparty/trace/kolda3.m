% KOLDA3 - Comptues [[C; U1, U2, U3]]
%
% Syntax
%  function X = kolda3(C, U1, U2, U3)
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


function X = kolda3(C, U1, U2, U3)

sz=[size(U1,1), size(U2,1), size(U3,1)];

dims=size(C);
if ndims(C)<3
  dims(3)=1;
end

U1=U1(:,1:dims(1));
U2=U2(:,1:dims(2));
U3=U3(:,1:dims(3));

for ii=1:sz(3)
  Cii=reshape(reshape(C,[dims(1)*dims(2), dims(3)])*U3(ii,:)',dims(1:2));
  X(:,:,ii)=U1*Cii*U2';
end
