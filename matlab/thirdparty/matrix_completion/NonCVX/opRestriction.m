function op = opRestriction(n,idx)
% OPRESTRICTION Restriction operator
%
%    OPRESTRICTION(N,IDX) creates a restriction operator that
%    selects the entries listed in IDX from an input vector of
%    length N. The adjoint of the operator creates a vector of
%    length N and fills the entries given by IDX with the input
%    data. 
%
%    See also opColumnRestriction, opMask.

%   Copyright 2007, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opRestriction.m 625 2007-10-23 02:50:10Z mpf $

if (min(idx) < 1) || (max(idx) > n)
   error('Index parameter must be integer and match dimensions of the operator');
end

op = @(x,mode) opRestriction_intrnl(n,idx,x,mode);

function y = opRestriction_intrnl(n,idx,x,mode)
if mode == 0
   m = length(idx);
   y = {m,n,[0,1,0,1],{'Restriction'}};
elseif mode == 1
   y = x(idx);
else
   y = zeros(n,1);
   y(idx) = x;
end
