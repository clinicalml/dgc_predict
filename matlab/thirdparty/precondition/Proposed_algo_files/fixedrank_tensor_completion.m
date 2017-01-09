function [Xsol, infos, problem_description] = fixedrank_tensor_completion(problem_description)
% Perform fixed-rank tensor completion with the proposed Riemannian metric.
%
% function Xsol = fixedrank_tensor_completion(problem_description)
% function [Xsol, costopt, infos, problem_description] = fixedrank_tensor_completion(problem_description)
%
% Input:
% -------
%
% problem_description: The problem structure with the description of the problem.
%
%
% - problem_description.data_train: Data structure for known entries that are used to learn a low-rank model.
%                                   It contains the 2 fields that are shown
%                                   below. An empty "data_train" structure
%                                   will generate a random instance.
%
%       -- data_train.entries:      A column vector consisting of known
%                                   distances. An empty "data_train.entries"
%                                   field will generate a random instance.
%
%
%       -- data_train.subs:         The [i j k] subscripts of corresponding
%                                   entries. An empty "data_train.subs"
%                                   field will generate a random instance.%
%
%
% - problem_description.data_test:  Data structure to the "unknown" (to the algorithm) entries.
%                                   It contains the 2 fields that are shown
%                                   below. An empty "data_test" structure
%                                   will not compute the test error.
%
%       -- data_test.entries:       A column vector consisting of "unknown" (to the algorithm)
%                                   entries. An empty "data_test.entries"
%                                   field will not compute the test error.
%       -- data_test.subs:          The [i j k] subscripts of the corresponding
%                                   entries. An empty "data_test.subs"
%                                   field will not compute the test error.
%
%
% - problem_description.data_validate:  Similar to "datat_test" description.
%                                       Any additional functionality.
%
%
% - problem_description.tensor_size: The size of the tensor. An empty
%                                    "tensor_size", but complete "data_train" structure
%                                    will lead to an error, to avoid
%                                    potential data inconsistency.
%
%
% - problem_description.tensor_rank: Rank. By default, it is [3 4 5].
%
%
%
% - problem_description.params:  Structure array containing algorithm
%                                parameters for stopping criteria.
%
%       -- params.tolgradnorm:   Tolerance on the norm of the gradient.
%                                By default, it is 1e-5.
%       -- params.tolmse:        Tolerance on the mlse.
%                                By default, it is 1e-5.
%       -- params.maxiter:       Maximum number of fixe-rank iterations.
%
%       -- params.solver:        Fixed-rank algorithm. Options are
%                                '@trustregions' for trust-regions,
%                                '@conjugategradient' for conjugate gradients,
%                                '@steepestdescent' for steepest descent.
%                                 By default, it is '@conjugategradient'.
%
%       -- params.linesearch:    Stepsize computation. Options are
%                                '@linesearchguess2' for guessed stepsize
%                                 with degree 2 polynomial guess,
%                                '@linesearchdefault' by Manopt
%                                 By default, it is '@linesearchguess2'.
%
%      -- params.show_plots:     Show basic plots. 
%                                By default, it is false; 
%
%      -- params.verbosity:      Show output with 0, 1, 2
%                                By default, 2; 
%
%
% Output:
% --------
%
%   X_sol:                Solution.
%   infos:                Structure array with computed statistics.
%   problem_description:  Structure array with used problem description.
%
%
% Please cite the Manopt paper as well as the research paper:
%     @TechReport{kasai2015precon,
%       Title   = {{R}iemannian preconditioning for tensor completion},
%       Author  = {Kasai, H. and Mishra, B.},
%       Journal = {Arxiv preprint arXiv:1506.02159},
%       Year    = {2015}
%     }


% Original authors: Hiroyuki Kasai and Bamdev Mishra, June 05, 2015.
% Contributors: 
% Change log:


        
    % Check problem description
    if ~exist('problem_description', 'var')
        problem_description = struct();
    end
    problem_description = check_problem_description(problem_description); % Check the problem description;
    
    
    % Define the problem data
    data_train = problem_description.data_train;
    data_train.nentries = length(data_train.entries); % Useful
    
    data_test =  problem_description.data_test; % Test data
    data_validate =  problem_description.data_validate; % Validation data
    
    tensor_size = problem_description.tensor_size;
    tensor_rank = problem_description.tensor_rank;
    params =  problem_description.params;
    
    n1 = tensor_size(1);
    n2 = tensor_size(2);
    n3 = tensor_size(3);
    r1 = tensor_rank(1);
    r2 = tensor_rank(2);
    r3 = tensor_rank(3);
    
    
    % The fixed-rank Tucker factory
    % Create the problem structure
    % quotient geometry
    M = fixedrankfactory_tucker_preconditioned(tensor_size, tensor_rank);
    problem.M = M;
    
    
    % Problem description: cost
    
    problem.cost = @cost;
    function [f store] = cost(X, store)
        if ~isfield(store, 'residual_vec') % Re use if it already computed
            store.residual_vec = compute_residual(X, data_train);
            residual_vec = store.residual_vec;
            store.cost = 0.5*(residual_vec'*residual_vec); % BM
        end
        f = store.cost;
    end
    
    
    
    % Problem description: gradient
    % We need to give only the Euclidean gradient, Manopt converts it to 
    % to the Riemannian gradient internally.
    problem.egrad =  @egrad; 
    function [g store] = egrad(X, store)
        if ~isfield(store, 'residual_vec') % Re use if it already computed
            store.residual_vec = compute_residual(X, data_train);
        end
        
       
        [temp,temp1,temp2,temp3] = calcProjection_mex(data_train.subs', store.residual_vec, X.U1', X.U2', X.U3' );
        
        
        Y23 = reshape(temp1, [n1 r2 r3]); % tensor(reshape(temp1, [n1 r2 r3]));
        Y13 = reshape(temp2, [r1 n2 r3]); % tensor(reshape(temp2, [r1 n2 r3]));
        Y12 = reshape(temp3, [r1 r2 n3]); % tensor(reshape(temp3, [r1 r2 n3]));
        
        T1 = reshape(Y23, n1, r2*r3) * reshape(X.G, r1, r2*r3)';
        T2 = reshape(permute(Y13, [2 1 3]), n2, r1*r3) * reshape(permute(X.G, [2 1 3]), r2, r1*r3)';
        T3 = reshape(permute(Y12, [3 1 2]), n3, r1*r2)* reshape(permute(X.G, [3 1 2]), r3, r1*r2)';

        g.U1 = T1;
        g.U2 = T2;
        g.U3 = T3;
        g.G = reshape(temp, [r1 r2 r3]);% tensor(reshape(temp, [r1 r2 r3]), [r1, r2, r3]);
    end
    
    
    % Problem description: stepsize computation with degree 2 polynomial
    % approximation
    function [t store] = linesearchguess2(X, eta, store)
        t  = stepsize_guess_degree2(X, eta, store); % Better for large-scale instances.
    end
    
    function tmin  = stepsize_guess_degree2(X, eta, store)
        if ~isfield(store, 'residual_vec') % Re use if it already computed
            store.residual_vec = compute_residual(X, data_train);
        end
        residual_vec = store.residual_vec;
        
        % new compute
        tmin = compute_stepsize_initial(X, eta, residual_vec, data_train.subs);
    end
    
    
    
    %     % Check numerically whether gradient and Hessian are correct
    %     checkgradient(problem); %# ok
    %     drawnow;
    %     pause;
    
    
    
    % Options
    
    
    % Ask Manopt to compute the error at every iteration when a
    % test and/or validate datasets are provided.
    if ~isempty(data_test) || ~isempty(data_validate)
        options.statsfun = @compute_test_validate_error;
    end
    function stats = compute_test_validate_error(problem, X, stats)
        if ~isempty(data_test)
            stats.test_error_square = compute_cost(X, data_test);
        end
        if ~isempty(data_validate) % If there is data_validate, then compute it at every iteration
            stats.validate_error_square = compute_cost(X, data_validate);
        end
    end
    
    
    % Call appropriate algorithm
    solver = params.solver; % Any algorithm that Manopt supports
    options.maxiter = params.maxiter;    
    options.tolgradnorm =  params.tolgradnorm;
    options.verbosity =  params.verbosity;
    
    
    % MSE stopping criteria options
    options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, Y, info, last) %#ok<INUSL>
        stopnow = (last >= 3 && 2*(info(last).cost)/data_train.nentries < params.tolmse);
    end
        
    
    if strcmp(func2str(params.linesearch), 'linesearchguess2')
        problem.linesearch = @linesearchguess2;
    elseif strcmp(func2str(params.linesearch), 'linesearchdefault')
        options.linesearch = @linesearch;
    else
        warning('fixedrank_tensor_completion:linesearch', ...
            'Linesearch is not properly defined. We work with default Manopt option.\n3');
        options.linesearch = @linesearch; 
    end
    
    % Call the solver
    X_initial = problem_description.X_initial;
    [Xsol, unused, infos] = solver(problem, X_initial, options);
    
    % Plots
    show_plots(problem_description, infos);
    
    
    
    % Third party Mex files by Michael Steinlechner <michael.steinlechner@epfl.ch> (MS)
    
    % Use of MS Mex file on computing entries of a tensor stored in Tucker format
    function vals = compute_residual(X, A_Omega)
        %CALCGRADIENT Calculate the euclid. gradient of the obj. function
        %   Wrapper function for calcGradient_mex.c
        %
        %   Computes the euclid. gradient of the objective function
        %
        %         X_Omega - A_Omega
        %
        %   between a sparse tensor A_Omega and a Tucker tensor X.
        
        %   GeomCG Tensor Completion. Copyright 2013 by
        %   Michael Steinlechner
        %   Questions and contact: michael.steinlechner@epfl.ch
        %   BSD 2-clause license, see LICENSE.txt
        
        vals = calcGradient_mex(A_Omega.subs', A_Omega.entries, ...
            X.G, X.U1', X.U2', X.U3');
        
    end
    
    function cost = compute_cost(X, A_Omega)
        %CALCFUNCTION Calculate the value of the objective function
        %   Wrapper function for calcFunction_mex.c
        %
        %   Computes the value of the objective Function
        %
        %       0.5 * || X_Omega - A_Omega ||^2
        %
        %   See also calcGradient
        
        %   GeomCG Tensor Completion. Copyright 2013 by
        %   Michael Steinlechner
        %   Questions and contact: michael.steinlechner@epfl.ch
        %   BSD 2-clause license, see LICENSE.txt
        
        cost = 0.5 * calcFunction_mex(A_Omega.subs', A_Omega.entries, ...
            X.G, X.U1', X.U2', X.U3');
    end
    
    
    
    function t = compute_stepsize_initial(X, eta, residual_vec, subs)
        %CALCINITIAL Calculate the initial guess for the line search.
        %
        %   Wrapper function for calcInitial_mex.c
        %
        %   See also calcGradient, calcFunction
        %
        
        %   GeomCG Tensor Completion. Copyright 2013 by
        %   Michael Steinlechner
        %   Questions and contact: michael.steinlechner@epfl.ch
        %   BSD 2-clause license, see LICENSE.txt
        
        t = -calcInitial_mex(subs', residual_vec, ...
            X.G, X.U1', X.U2', X.U3',...
            eta.G, eta.U1', eta.U2', eta.U3'); 
    end

    
    
end






%% Local defaults
function localdefaults = getlocaldefaults()
    localdefaults.tolgradnorm = 1e-5;
    localdefaults.tolmse =1e-5;
    localdefaults.maxiter = 250;
    localdefaults.solver = @conjugategradient; % Conjugate gradients
    localdefaults.linesearch = @linesearchguess2; % Stepsize guess with 2 polynomial
    localdefaults.show_plots = false;
    localdefaults.verbosity = 2;
end




%% Problem description check
function checked_problem_description = check_problem_description(problem_description)
    
    checked_problem_description = problem_description;

    
    % Check train data
    if isempty(problem_description)...
            || ~all(isfield(problem_description,{'data_train'}) == 1)...
            || ~all(isfield(problem_description.data_train,{'entries', 'subs'}) == 1)...
            || isempty(problem_description.data_train.entries)...
            || isempty(problem_description.data_train.subs)
        
        warning('fixedrank_tensor_completion:problem_description', ...
            'The training set is either empty or not properly defined. We work with a random instance.\n');
        checked_problem_description = get_random_instance();
        return; % No need for further check
    end
    
    
    % Check tensor size
    if ~isfield(problem_description, 'tensor_size')
        error('fixedrank_tensor_completion:problem_description',...
            'Error. The dimensions of the tensor corresponding to field "tensor_size" must be given. \n');
    end
    
    
    % Check tensor rank
    if ~isfield(problem_description, 'tensor_rank')...
            || isempty(problem_description.tensor_rank)
        warning('fixedrank_tensor_completion:problem_description', ...
            'The field "tensor_rank" is properly defined. We work with "[3 4 5]".\n');
        tensor_rank = [3 4 5];
    else
        tensor_rank = problem_description.tensor_rank;
    end
    checked_problem_description.tensor_rank = tensor_rank;
   
    
    
    % Check initialization
    if ~isfield(problem_description, 'X_initial')...
            || isempty(problem_description.X_initial)
        warning('fixedrank_tensor_completion:problem_description', ...
            'The field "X_initial" is not given. We work with a random initialization.\n');
        X_initial = [];
    else
        X_initial = problem_description.X_initial;
    end
    checked_problem_description.X_initial = X_initial;
    

    
    % Check testing dataset
    if ~isfield(problem_description,{'data_test'})...
            || ~all(isfield(problem_description.data_test,{'entries', 'subs'}) == 1)...
            || isempty(problem_description.data_test.entries)...
            || isempty(problem_description.data_test.subs)
        
        warning('fixedrank_tensor_completion:problem_description', ...
            'The field "data_test" is either empty or not properly defined. We work with the default "[]".\n');
        data_test = [];
    else
        data_test = problem_description.data_test;
    end
    checked_problem_description.data_test = data_test;
    
    
    
    % Check validation dataset, an extra set that may be used
    if ~isfield(problem_description,{'data_validate'})...
            || ~all(isfield(problem_description.data_validate,{'entries', 'subs'}) == 1)...
            || isempty(problem_description.data_validate.entries)...
            || isempty(problem_description.data_validate.subs)
        
        warning('fixedrank_tensor_completion:problem_description', ...
            'The field "data_validate" is either empty or not properly defined. We work with the default "[]".\n');
        data_validate = [];
    else
        data_validate = problem_description.data_validate;
    end
    checked_problem_description.data_validate = data_validate;
    
    
    % Check parameters
    if isfield(problem_description, 'params')
        params = problem_description.params;
    else
        params = struct();
    end
    params = mergeOptions(getlocaldefaults(), params);
    checked_problem_description.params = params;
end


%% A radom problem instance
function problem_description = get_random_instance()
    
    % tensor size
    tensor_size = [300 200 100];
    
    % core tensor dimensions
    core_dims = [3, 4, 5];
    
    
    % Sampling ratio
    OS = 10;
    dim = tensor_size(1)*core_dims(1) - core_dims(1)^2 ...
        + tensor_size(2)*core_dims(2) - core_dims(2)^2 ...
        + tensor_size(3)*core_dims(3) - core_dims(3)^2 ...
        + core_dims(1)*core_dims(2)*core_dims(3);
    
    fraction = (OS*dim)/(tensor_size(1)*tensor_size(2)*tensor_size(3));
    
    % Noise level
    noise = [0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1];
    noise_level = 1;
    
    
    % A random tucker tensor
    Xttensor = makeRandTensor(tensor_size, core_dims );
    
    
    % Generate a training set
    subs = makeOmegaSet(tensor_size, round(fraction * prod(tensor_size)));% Training indices sorted as suggested by MS
    vals = getValsAtIndex_mex(subs', Xttensor.G, Xttensor.U1', Xttensor.U2', Xttensor.U3'); % Training entries
    nentries = size(subs, 1); % Training number of known entries
    noise_vector = rand(nentries, 1);
    inverse_snr = noise(noise_level);
    vals_noisy = vals + (inverse_snr * norm(vals) / norm(noise_vector))*noise_vector; % entries added with noise
    
    data_train.entries = vals_noisy;
    data_train.subs = subs;
    data_train.nentries = nentries;
    
    
    % Generate a test set
    data_test.nentries = 1*data_train.nentries; % Depends how big data can we handle
    data_test.subs = [randi(tensor_size(1), data_test.nentries, 1), ...
        randi(tensor_size(2), data_test.nentries, 1), ...
        randi(tensor_size(3), data_test.nentries, 1)];
    
    data_test.entries = getValsAtIndex_mex(data_test.subs', Xttensor.G, Xttensor.U1', Xttensor.U2', Xttensor.U3');
    
    
    % Generate a validate set, an extra set that may be used
    data_validate = [];
    
    
    % Rank
    tensor_rank = core_dims;
    
    
    
    % Basic parameters used in optimization
    params = struct();
    params = mergeOptions(getlocaldefaults, params);
    params.show_plots = true;
    
    
    % Problem description
    problem_description.data_train = data_train;
    problem_description.data_test = data_test;
    problem_description.data_validate = data_validate;
    problem_description.tensor_size = tensor_size;
    problem_description.tensor_rank = tensor_rank;
    problem_description.X_initial = [];
    problem_description.params = params;
    
end



function show_plots(problem_description, infos)
    
    if ~problem_description.params.show_plots;
        return;
    end
    
    
    solver = problem_description.params.solver;
    linesearch = problem_description.params.linesearch;
    
    data_test = problem_description.data_test;
    nentries = length(problem_description.data_test.entries);
    fraction = nentries/prod(problem_description.tensor_size);
    
    % Training on Omega
    fs = 20;
    figure;
    semilogy([infos.iter], [infos.cost], '-O','Color','blue','linewidth', 2.0);
    ax1 = gca;
    set(ax1,'FontSize',fs);
    xlabel(ax1,'Iterations','FontSize',fs);
    ylabel(ax1,'Cost','FontSize',fs);
    legend([func2str(solver),' + ', func2str(linesearch)]);
    legend 'boxoff';
    box off;
    title(['Fraction known ',num2str(fraction),', Training'])
    
    if ~isempty(data_test)
        % Testing
        fs = 20;
        figure;
        semilogy([infos.test_error_square], '-O','Color','blue','linewidth', 2.0);
        ax1 = gca;
        set(ax1,'FontSize',fs);
        xlabel(ax1,'Iterations','FontSize',fs);
        ylabel(ax1,'Test error','FontSize',fs);
        legend([func2str(solver),' + ',func2str(linesearch)]);
        legend 'boxoff';
        box off;
        title(['Fraction known ',num2str(fraction),', test error on a set  \Gamma, different from \Omega'])
        
    end 
end


