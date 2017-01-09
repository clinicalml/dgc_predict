function [M] = Unfold( T, dim, i )
M = reshape(shiftdim(T,i-1), dim(i), []);