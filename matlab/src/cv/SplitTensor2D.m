function [testIdx, trainIdx, nIter] = SplitTensor2D(T, ...
    nsplit, maxsplit, minSigsPerDrug, printFlag, maxTries)
% DESCRIPTION: This function randomly splits data into nsplit folds (if the
% data does not divide evenly, some signatures will not make it into a
% fold)
%
% INPUTS:
% T         Tensor
% nsplit    Number of folds to compute
% maxsplit  Number of these folds (<= nsplit) to actually output
% minSigsPerDrug 
% printFlag 
% 
% OUTPUTS: 

% testIdx   Cell array of length maxsplit containing linear indices into
%           drug-cell matrix specifying test data
% trainIdx  Cell array of length maxsplit containing linear indices into
%           drug-cell matrix specifying training data
% nIter     number of random samples before the folds were balanced



%% set defaults
if ~exist('nsplit') || length(nsplit) == 0
    nsplit = 10;
end

if ~exist('maxsplit') || length(maxsplit) == 0
    maxsplit = nsplit;
end

if ~exist('minSigsPerDrug')
    minSigsPerDrug = 1;
end

if ~exist('printFlag')
    printFlag = false;
end

if ~exist('maxTries')
  maxTries = 20000;
end

maxsplit = min(maxsplit, nsplit);

%% find locations of all drug/cell combinations that contain data
A = squeeze(~isnan(T(:,1,:)));
I_obs = find(A);


%% draw a random permutation corresponding to this set of drug/cell indices
N = length(I_obs);
sz = size(T);


%% finally, define folds

% round down to divide nsplit evenly
n = N - mod(N,nsplit); 

% number of signatures per fold
k = n / nsplit;

trainIdx = cell(maxsplit, 1);
testIdx = cell(maxsplit, 1);

balanced = false;
nIter = 1;

while ~balanced && nIter < maxTries
    p = randperm(N);

    for i = 1:maxsplit
        idx1 = (i-1)*k + 1;
        idx2 = i*k;
        testIdx{i} = I_obs(p(idx1:idx2));
        trainIdx{i} = setdiff(I_obs, testIdx{i});
    end
    
    balanced = IsBalanced(trainIdx, minSigsPerDrug, sz);
    
    if ~balanced
        if printFlag
            disp('Split is unbalanced, reshuffling...')
        end
        nIter = nIter + 1;
    end
    
    if nIter == maxTries
        error(sprintf('Could not find a balanced split in %s tries.', maxTries));
    end
    if mod(nIter,1000)==0
        fprintf('%d attempts at tensor split\n', nIter)
    end
end

end






