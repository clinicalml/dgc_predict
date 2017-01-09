function d = ComputeDensity(T)
% Computes the observation density

d = length(find(~isnan(T))) / prod(size(T));
