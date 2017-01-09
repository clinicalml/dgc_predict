
% arguments to pftest
options(1) = 1e-5; % convergence criterion
options(2) = 2; % random initialization of missing values
options(3) = 0; % amount of plots
options(4) = 2; % no scaling (vs options(4)=1: default scaling- columns in mode one carry the variance)
options(5) = 10; % how often to show fit
options(6) = 1000; % max number of iterations
numReps = 3;
maxFactors = 4;

baseDir = '/Users/rhodos/Desktop/Dropbox/LINCS/data/expr/tensor/';
files = dir(baseDir);
idx = 6:3:15;

err_nonorm = cell(length(idx));
cor_nonorm = cell(length(idx));
iter_nonorm = cell(length(idx));


err_norm = cell(length(idx));
cor_norm = cell(length(idx));
iter_norm = cell(length(idx));

c = 1;
for i = 10 %4:length(files) % first three are '.', '..', and 'README'
   disp(files(i).name)
   load([baseDir files(i).name]) 
   
   disp(sprintf('num perts = %d', size(T,1)));
   
   %[err_norm{c},cor_norm{c},iter_norm{c}] = pftest(numReps,T,maxFactors,options);
   
   % "un" normalize data
   T = T .* repmat(sqrt(v), size(T,1), 1, size(T,3));
   T = T + repmat(m, size(T,1), 1, size(T,3));
   [err_nonorm{c},cor_nonorm{c},iter_nonorm{c}] = pftest(numReps,T,maxFactors,options);
   
   c = c+1;
end




