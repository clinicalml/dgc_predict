function args = GetArgs(model, T_testIdx, testIdx, sz)

if ~exist('T_testIdx')
    T_testIdx = [];
end

if ~exist('testIdx')
    testIdx = [];
end

if ~exist('sz')
    sz = [];
end

switch(model)
    case 'tmac'
        args = GetTmacArgs(T_testIdx, testIdx);
    case 'asmatrix'
        args = GetAsMatrixArgs();
    case 'constrained'
        args = GetConstrainedArgs();
    case 'mixture'
        args = GetMixtureArgs();
    case 'ha_lrtc'
        args = GetHaLRTCArgs(T_testIdx, testIdx);
    case 'si_lrtc'
        args = GetSiLRTCArgs(T_testIdx, testIdx);
    case 'fa_lrtc'
        args = GetFaLRTCArgs(sz, T_testIdx, testIdx);
    case 'mean' % 1D-Mean
        args = [];
    case 'mean2'% 2D-Mean
        args = GetMean2DArgs();
    case 'matrixcomp'
        args = GetMatrixCompArgs();
    case 'knnd'
        args = GetKNNDrugArgs();
    case 'knnc'
        args = GetKNNCellArgs();
    case 'knndc'
        args = GetKNNDrugCellArgs();
    case 'unfoldmc'
        args = GetUnfoldMCArgs();
    otherwise
        error('invalid model');
end
end

function args = GetTmacArgs(T_testIdx, testIdx)

args.estCoreNway = [10 10 3];

% setup options
opts.maxit = 100;
opts.tol = 1e-4;
opts.alpha = [1 1 1];
opts.alpha_adj = 0;
opts.rank_adj = [1 1 1]; 
opts.rank_max = [50 30 6]; 
opts.rank_inc = [1 1 1];   
opts.rank_min = [1 1 1];
opts.T_testIdx = T_testIdx;
opts.testIdx = testIdx;

args.opts = opts;

end

function args = GetAsMatrixArgs() 
alpha = [1 0 1];
args.alpha = alpha / sum(alpha);
args.eta = 32;
args.tol = 5e-2;
end

function args = GetConstrainedArgs()
args.lambda = 0;
args.eta = 32;
args.tol = 5e-2;
end

function args = GetMixtureArgs()
args.lambda = 0;
args.eta = 0.1;
args.tol = 1e-5;
end

function args = GetHaLRTCArgs(T_testIdx, testIdx)
args.Ttest = T_testIdx;
args.Omega_test = testIdx;
alpha = [1 1e-3 1];
args.alpha = alpha / sum(alpha);
args.beta = 1;
args.maxIter = 500;
args.epsilon = 5e-2;
end

function args = GetSiLRTCArgs(T_testIdx, testIdx)
args.Ttest = T_testIdx;
args.Omega_test = testIdx;
alpha = [1 1e-3 1];
args.alpha = alpha / sum(alpha);
args.beta = 32*[1 1 1];
args.maxIter = 500;
args.epsilon = 1e-10;
end

function args = GetFaLRTCArgs(sz, T_testIdx, testIdx)
args.Mtest = T_testIdx;
args.Omega_test = testIdx;
alpha = [1 1e-2 1];
args.alpha = alpha / sum(alpha);
args.mu = 0.01 * alpha ./ sqrt(sz);
args.C =  0.5;
args.L0 = 1e-5;
args.maxIter = 20;
args.epsilon = 1e-1; 
end

function args = GetMean2DArgs()
args.alpha = 0.5;
end

function args = GetMatrixCompArgs()
args.p = 1;
end

function args = GetKNNDrugArgs()
args.K = 10;
% if exist('pair') && length(pair) == 2
%     args.d = pair(1);
%     args.c = pair(2);
% end
end

function args = GetKNNCellArgs()
args.K = 10;
end

function args = GetKNNDrugCellArgs()
args.knnd = GetKNNDrugArgs();
args.knnc = GetKNNCellArgs();
%args.alpha = 0.5;
end

function args = GetUnfoldMCArgs()
args.p = 1;
args.unfold_dim = 1;
end
