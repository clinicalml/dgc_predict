function k = NumSigs(T, dim)
% Computes number of signatures in the tensor. Assumes that if the
% expression signature is available, it will be either completely measured,
% or completely absent.

A = squeeze(T(:,1,:));


if ~exist('dim')
    k = sum(~isnan(A(:)));
elseif strcmp(dim, 'drug')
    k = sum(~isnan(A),2);
elseif strcmp(dim, 'cell')
    k = sum(~isnan(A),1);
else
    error('unexpected value for dim');
end

end