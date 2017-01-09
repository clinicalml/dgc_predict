function T = RandTensor(dims, k, noise)
% generate random tensor of specified rank k, following parafac model
% density defines the observation density, which will be met be removing
% entire columns from the second dimension (i.e. gene signatures).

if ~exist('dims')
    dims = [20 20 20];
end

if ~exist('k')
    k = 3;
end

if ~exist('noise')
    noise = 0;
end

T = zeros(dims);

for kk = 1:k
    x = randn(dims(1), 1);
    y = randn(dims(2), 1);
    z = randn(dims(3), 1);
    T = T + OuterProduct3(x,y,z);
end

if noise > 0
    T = T + randn(dims(1), dims(2), dims(3))*noise; 
end

end