function [T_model, out] = CompleteTensorKNNDrugCell(T, args)

if ~isfield(args, 'alpha')     % run inner loop to choose alpha, just holding out 20% of the data
    disp('Choosing alpha in inner loop..');
    
    % subset to cases where I have at least 3 signatures per drug, so as to
    % avoid issues with splitting
    Tdense = T(NumSigs(T, 'drug') >= 3,:,:);
    assert(size(Tdense, 1) > 0.5*size(T,1)); 
    
    % first get estimates from KNNc and KNNd
    T_est = TensorCV4({'knnc','knnd'}, Tdense, 'test', 5, 1);
    
    % then try various linear combinations
    alpha_vec = 0:0.05:1;
    PCT = nan(length(alpha_vec),1);
    for a = 1:length(alpha_vec)
        alpha = alpha_vec(a);
        T_tmp = T_est{1}*alpha + T_est{2}*(1-alpha);
        PCT(a) = corr(T_tmp(:), Tdense(:), 'rows', 'pairwise');
    end
    
    % then find the best combination
    [~, best_idx] = max(PCT);
    out.alpha_best = alpha_vec(best_idx);
    out.PCT = PCT;
    alpha = out.alpha_best;
    fprintf('..setting alpha=%0.2f\n', alpha);

else % or if alpha is defined in args, just use that

    disp('Using input alpha');
    alpha = args.alpha;
    out = [];

end

T1 = CompleteTensorKNNCell(T, args.knnc);
T2 = CompleteTensorKNNDrug(T, args.knnd);
T_model = alpha * T1 + (1-alpha)*T2;

end
