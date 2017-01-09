
InitRand(123);

models = {'knnd', 'knnc', 'fa_lrtc'};
T = RandTensor([1,1,1]*20,4,0.1);
exp_name = 'test';
nFolds = 10;
maxFolds = 10;
saveFile = false;
printFlag = false;
debugFlag = false;
normalize = true;

[T_imputed, PCT, PCTf, time, testIdx] = TensorCV4(models, T, ...
    exp_name, nFolds, maxFolds, saveFile, printFlag, debugFlag, ...
    normalize);

assert(mean(PCT) > 0.5);

disp('Warning: need to do further testing on TensorCV4');
% look at PCT
% mean(PCT, 1)
% Hm not sure exactly why mean and mean2 are getting negative correlations.
% I think it's because the model mismatch is too strong.

% look at locations/numbers of missing elements in completed tensor
%T_missing = RemoveSigsRandom(T, 100, true);    
%[T_imputed, PCT, time, testIdx] = TensorCV4(models, T_missing);

