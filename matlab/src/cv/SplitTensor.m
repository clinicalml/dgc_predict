function [testIdx, trainIdx] = SplitTensor(T, nsplit, maxsplit, pertIdx)
% DESCRIPTION: This function randomly splits data into nsplit folds (if the
% data does not divide evenly, some signatures will not make it into a
% fold)
%
% INPUTS:
% T         Tensor
% nsplit    Number of folds to compute
% maxsplit  Number of these folds (<= nsplit) to actually output
% pertIdx   Optional argument used to subset tensor
% 
% OUTPUTS:
% testIdx   Cell array of length maxsplit, where each entry contains linear
%           indices into T specifying test data
% trainIdx  Cell array of length maxsplit, where each entry contains linear
%           indices into T specifying training data


%% set defaults
if ~exist('nsplit') || length(nsplit) == 0
    nsplit = 10;
end

if ~exist('maxsplit') || length(maxsplit) == 0
    maxsplit = nsplit;
end
        
maxsplit = min(maxsplit, nsplit);

%% find locations of all drug/cell combinations that contain data
A = squeeze(~isnan(T(:,1,:)));
I_full = find(A);

%% subset using pertIdx if it is defined
if exist('pertIdx')
    nanIdx = setdiff(1:size(A,1), pertIdx);
    A(nanIdx,:) = 0;
    I = find(A);
else
   I = I_full; 
end

%% draw a random permutation corresponding to this set of drug/cell indices
n = length(I);
p = randperm(n);


%% finally, define folds

% round down to divide nsplit evenly
n = n - mod(n,nsplit); 

% number of signatures per fold
k = n / nsplit;

testIdx = cell(maxsplit, 1);
trainIdx = cell(maxsplit, 1);

for i = 1:maxsplit
    idx1 = (i-1)*k + 1;
    idx2 = i*k;
    
    idxTest_linear2d = I(p(idx1:idx2));
    idxTrain_linear2d = setdiff(I_full, idxTest_linear2d);
    
    testIdx{i} = MapMatrixInd2Tensor(idxTest_linear2d,size(A), size(T));
    trainIdx{i} = MapMatrixInd2Tensor(idxTrain_linear2d, size(A), size(T));
end

end