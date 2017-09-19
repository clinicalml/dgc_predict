
% knnd_ts is the KNNDrug method (i.e. DNPP) where instead of using all cell
% types to compute the similarity matrix between pairs of drugs, one uses
% only the Top L Cell types with the highest similarity. This was based on
% a PSB reviewer's suggestion. However, we find that the results are a bit
% worse than just using the vanilla KNNDrug.

models = {'knnd_ts'};
exp_name = 'test';
nFolds = 10;
maxFolds = 10;
saveFile = false;
printFlag = true;
debugFlag = false;
normalize = true;
minSigsPerDrug = 1;

T = GetTensor('../../../data/expr/tensor/tsize/small/small');
 
maxPCTf = -Inf;
Lvec = 1:15;

PCTf = cell(length(Lvec), 1);
PCT = zeros(length(Lvec), 1);
time = zeros(length(Lvec), maxFolds);
args.K = 10;

for L = Lvec
    
    disp(L);
    
    args.maxSim = L;
    
    [T_est, PCT(L), PCTf{L}, time(L,:)] = TensorCV('knnd_ts', T,...
        exp_name, nFolds, maxFolds, saveFile, printFlag, ...,
        debugFlag, normalize, args, minSigsPerDrug);
    
    if mean(PCTf{L}) > maxPCTf
        maxPCTf = mean(PCTf{L});
        T_keep = T_est;
        PrintIf(sprintf('new best PCT (mean across %d folds) = %0.2f', ...
            nFolds, maxPCTf), printFlag);
    end
end

% PCT =
% 
%     0.6291
%     0.6361
%     0.6398
%     0.6410
%     0.6418
%     0.6422
%     0.6426
%     0.6430
%     0.6428
%     0.6427
%     0.6426
%     0.6426
%     0.6426
%     0.6426
%     0.6426

% vs. PCT for knnd: 0.643

