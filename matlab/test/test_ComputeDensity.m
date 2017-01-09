

T = randn(3,4,5);

T(1,:,:) = NaN;

assert(ComputeDensity(T) == (2/3));