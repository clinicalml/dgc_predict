function [T_imputed, PCT, PCTf, time, testIdx] = TensorCV4(models, T, ...
    exp_name, nFolds, maxFolds, saveFile, printFlag, debugFlag, ...
    normalize, args_in, minSigsPerDrug)

% DESCRIPTION: This function performs k-fold cross-validation of tensor
% completion using one or more tensor completion algorithms.
%
% INPUTS:
% models        List of tensor completion algorithm names (see GetArgs).
% T             Data tensor
% exp_name      name of experiment for directory creation
% nFolds        # of cross-validation folds.  Use -1 for leave-one-out.
% maxFolds      Maximum number of folds to compute, e.g. could divide data
%               into nFolds=10 but then only run calculation on maxFolds=2
%               of these.
% saveFile      Boolean; whether to save results to file.
% printFlag     Boolean; just prints some updates.
% debugFlag     Boolean; if true, just set missing tensor values to 1.
% normalize     Boolean; (set to true in paper). Normalizes each gene
%               expression profile to have unit length.
% args_in       Tensor completion arguments specific to the selected model.
%               This is passed to CompleteTensor. Normally 'args' is pulled
%               from the GetArgs function, but this option allows you to
%               input them yourself.
% minSigsPerDrug Minimum number of signatures per drug per training fold
% 
% OUTPUTS:
% T_imputed     Cross-validated predictions
% PCT           A vector of full-tensor PCT values computed for each model.
% PCTf          An array of PCT (Pearson correlation with truth) computed for each fold and model.
% time          Runtime (seconds) per fold per model.
% testIdx       See 'SplitTensor'.

%% initialize, set defaults and load data
InitRand();

if ischar(models)
    models = {models};
end

if ~exist('exp_name')
    exp_name = 'test';
end

if ~exist('nFolds')
    nFolds = 10; 
end

if nFolds == -1 % code for leave-one-out
    nFolds = NumSigs(T);
end

if ~exist('maxFolds')
    maxFolds = min(nFolds, 10);
end

if ~exist('saveFile')
    saveFile = false;
end

if ~exist('printFlag')
    printFlag = true;
end

if ~exist('debugFlag')
    debugFlag = false;
end

if ~exist('normalize')
    normalize = true;
end

if ~exist('minSigsPerDrug')
    minSigsPerDrug = 1;
end


%% split data into folds
PrintIf(sprintf('splitting tensor...'), printFlag);
[testIdx, trainIdx, nIter] = SplitTensor2D(T, nFolds, maxFolds, minSigsPerDrug);
sigsPerFold = length(testIdx{1});
PrintIf(sprintf('Holding out %d sigs per fold, running %d out of %d possible folds\n', ...
    sigsPerFold, maxFolds, nFolds), printFlag);

%% initialize outputs
nModels = length(models);
PCT = nan(1, nModels);
PCTf = nan(maxFolds, nModels);
time = nan(maxFolds, nModels);
T_imputed = cell(nModels, 1); 
meas_full = T(:);

%% setup output files
model_str = StrJoin(models, '__');
outDir = sprintf('%sresults/%s/%ddrugs_%dcells/', DataDir(), exp_name, size(T,1), size(T,3));
exp_desc = sprintf('%s_%dfoldCV', model_str, nFolds);

%% create output directories if they do not exist
if ~exist(outDir) && saveFile
    mkdir(outDir)
end

%% for each model, run maxFolds folds and compute accuracy on held out data
for m = 1:nModels
    PrintIf(sprintf('Running %s model', models{m}), printFlag);
    T_imputed{m} = nan(size(T));
    
    for f = 1:maxFolds        
        
        PrintIf(sprintf('..fold %d', f), printFlag);
       
        % construct 'training' tensor containing all data except one fold 
        T_train = ConstructTrainingTensorFrom2D(T, trainIdx{f});
        
        % Get 3D test indices for fold
        testIdx3D = MapMatrixInd2Tensor(testIdx{f}, [size(T,1), size(T,3)], size(T));
        
        % check that there are enough signatures per drug
        min_sigs_per_drug = min(NumSigs(T_train,'drug')); 
        if min_sigs_per_drug == 1
            disp('Warning: some drugs have only one signature in this fold');
        else
            assert(min_sigs_per_drug >= minSigsPerDrug);
        end
             
        % get model-specific arguments
        args = GetArgs(models{m}, T(testIdx3D), testIdx3D, size(T));
        if exist('args_in')
            for fn = fieldnames(args_in)'
                args.(fn{1}) = args_in.(fn{1});
            end
        end
        
        % run tensor completion
        [T_model, time(f,m)] = CompleteTensor(T_train, models{m}, args, ...
                                         printFlag, debugFlag, normalize);
        
        % save predictions in tensor
        T_imputed{m}(testIdx3D) = T_model(testIdx3D);

        % compute accuracy per fold (PCTf) 
        meas = T(testIdx3D);
        est = T_model(testIdx3D);
        PCTf(f,m) = corr(meas, est);
        
        clear args
    end
    
    % compute accuracy on full tensor results
    PCT(m) = corr(meas_full, T_imputed{m}(:), 'rows', 'pairwise'); 
    
    % print mean of errors across folds
    PrintIf(sprintf('RESULTS for %s: PCTf = %0.3f +/- %0.3f', ...
        models{m}, mean(PCTf(:,m)), std(PCTf(:,m))), printFlag);
end

% write to file
if saveFile
    save([outDir exp_desc], 'T_imputed', 'PCTf', 'PCT', 'time',...
        'testIdx', 'nIter');
end

end



