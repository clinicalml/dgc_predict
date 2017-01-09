function InitRand(seed)
% Initializes random number generator either with a fixed seed (so that
% results can be reproducible) or a user-specified seed.

if ~exist('seed')
    seed = 123;
end

rng('default');
rng(seed);

end