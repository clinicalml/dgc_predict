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
    case 'inftucker'
        args = GetInfTuckerArgs();
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

function args = GetInfTuckerArgs()
args.data = {};
args.data.beta = 2;

args.model = [];
args.model.dim = 3;
args.model.loss = 'Gauss';
args.model.em_total = 20;
args.model.var_total = 20;
args.model.verbose = 1;
args.model.process = 'GP';
args.model.lambda = 1;
args.model.sigma2 = 1e-5;

args.optimizer = [];
args.optimizer.algorithm = @L1General2_PSSgb;
args.optimizer.options.order = -1;
args.optimizer.options.maxIter = 20;
args.optimizer.options.verbose = 1;

args.uparam = [];
args.uparam.alpha = 1;
args.uparam.decision = 1e-8;
args.uparam.delta = 10;
args.uparam.distance = 'Euclidean';
args.uparam.gamma = 1;
args.uparam.kreg = 1; %1e-3?
args.uparam.norm_reg = 1;
args.uparam.nu = 0.5;
args.uparam.phi = 1;
args.uparam.savestorage = 0;
args.uparam.taper = 'Spherical';
args.uparam.tau = 1;
args.uparam.theta = 1;
args.uparam.tol = 1e-3;
args.uparam.type = 'rbf';
args.uparam.weight = 1;
end
