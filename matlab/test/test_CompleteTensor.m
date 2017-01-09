
% This tests how well these models recover entries of a rank-1 CP tensor.
% Other models (e.g. Mean-1D/2D, KNNs) do not model data as CP, so could
% not expect convergence.

InitRand(100)
TT = RandTensor([10 10 10], 1, 0.1);
d = [0.2 0.95]; 
ind_keep = VaryDensity(TT, d);
printFlag = true;

models = {'mixture', 'ha_lrtc', 'fa_lrtc', 'si_lrtc', 'constrained',...
    'asmatrix' 'tmac', 'matrixcomp', 'unfoldmc'};

err = zeros(length(d),1);

for m = 1:length(models)
    tic;
    model = models{m};
    PrintIf(sprintf('MODEL = %s', model), printFlag)
    for i = 1:length(d)
        args = GetArgs(model, [], [], size(TT));
        T = ConstructTrainingTensor(TT, ind_keep{i});
        
        if strcmp(model, 'knn')
            args.K = 1;
        end
        
        expression = 'T_out = CompleteTensor(T,model,args,printFlag);';
        evalc(expression);
        
        idx = find(isnan(T));
        err(i) = norm(T_out(idx) - TT(idx)) / sqrt(length(idx));
        
        PrintIf(sprintf('  density=%0.2f, err=%0.5f', d(i), err(i)), printFlag);    
    end
    assert(issorted(-err));
    time(m) = toc;
end

clear err

