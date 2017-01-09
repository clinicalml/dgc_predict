clc; clear all; close all;

%% Parameters

% tensor size
tensor_dims = [300 400 500]; %[3000 4000 5000];

% core tensor dimensions
core_dims = [12, 4, 3];


% Sampling ratio
OS = 5;
dim = tensor_dims(1)*core_dims(1) - core_dims(1)^2 ...
             + tensor_dims(2)*core_dims(2) - core_dims(2)^2 ...
             + tensor_dims(3)*core_dims(3) - core_dims(3)^2 ...
             + core_dims(1)*core_dims(2)*core_dims(3);

fraction = (OS*dim)/(tensor_dims(1)*tensor_dims(2)*tensor_dims(3));

% Noise level
noise = [0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1];
noise_level = 1;


% Tolerance
tolmse = 1e-8;
tolgradnorm = 1e-8;
maxiter = 200;



%% Generate data

Xttensor = makeRandTensor(tensor_dims, core_dims ); % Original tensor




%% Generate the training set

subs = makeOmegaSet(tensor_dims, round(fraction * prod(tensor_dims)));  % Training indices sorted as suggested by MS
vals = getValsAtIndex_mex(subs', Xttensor.G, Xttensor.U1', Xttensor.U2', Xttensor.U3'); % Training entries
nentries = size(subs, 1); % Training number of known entries
noise_vector = rand(nentries, 1);
inverse_snr = noise(noise_level);
vals_noisy = vals + (inverse_snr * norm(vals) / norm(noise_vector))*noise_vector; % entries added with noise

data_train.entries = vals_noisy; % BM
data_train.subs = subs; % BM
data_train.size = tensor_dims; % BM
data_train.nentries = nentries; % BM








%% Generate a test set

data_test.nentries = 1*data_train.nentries; % Depends how big data can we handle
data_test.subs = [randi(tensor_dims(1), data_test.nentries, 1), ...
    randi(tensor_dims(2), data_test.nentries, 1), ...
    randi(tensor_dims(3), data_test.nentries, 1)];

data_test.entries = getValsAtIndex_mex(data_test.subs', Xttensor.G, Xttensor.U1', Xttensor.U2', Xttensor.U3');
data_test.size = tensor_dims;


% If you dont want to compute an error on the test data, 
% then set data_test = [].
% data_test = [];



%% Generate a validate set

data_validate.nentries = 1*data_train.nentries; % Depends how big data can we handle
data_validate.subs = [randi(tensor_dims(1), data_validate.nentries, 1), ...
    randi(tensor_dims(2), data_validate.nentries, 1), ...
    randi(tensor_dims(3), data_validate.nentries, 1)];

data_validate.entries = getValsAtIndex_mex(data_validate.subs', Xttensor.G, Xttensor.U1', Xttensor.U2', Xttensor.U3');
data_validate.size = tensor_dims;


% If you dont want to compute an error on the validate data, 
% then set data_validate = [].
data_validate = [];




%%
X_init = makeRandTensor( tensor_dims, core_dims);



%% Call to our algorithm 


%  Required problem description
problem_description.data_train = data_train;
problem_description.data_test = data_test;
problem_description.data_validate = data_validate;
problem_description.tensor_size = tensor_dims;
problem_description.tensor_rank = core_dims;

problem_description.X_initial.U1 = X_init.U1;
problem_description.X_initial.U2 = X_init.U2;
problem_description.X_initial.U3 = X_init.U3;
problem_description.X_initial.G = X_init.G; 


% Some options, but not mandatory
problem_description.params.tolgradnorm = tolgradnorm; 
problem_description.params.tolmse = tolmse;
problem_description.params.maxiter = maxiter;
problem_description.params.solver = @conjugategradient; % Conjugate gradients
problem_description.params.linesearch = @linesearchguess2; % Stepsize guess with 2 polynomial


% Alogorithm
[Xsol, infos] = fixedrank_tensor_completion(problem_description);



%% Plots
solver = problem_description.params.solver;
linesearch = problem_description.params.linesearch;


% Training on Omega
fs = 20;
figure;
semilogy([infos.time], [infos.cost], '-O','Color','blue','linewidth', 2.0);
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel(ax1,'Time in seconds','FontSize',fs);
ylabel(ax1,'Train square error','FontSize',fs);
legend('Proposed algo');
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
    ylabel(ax1,'Test square error','FontSize',fs);
    legend('Proposed algo');
    legend 'boxoff';
    box off;
    title(['Fraction known ',num2str(fraction),', test error on a set  \Gamma, different from \Omega'])
    
end


