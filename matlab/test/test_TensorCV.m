
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

[T_imputed, PCT, PCTf, time, testIdx] = TensorCV(models, T, ...
    exp_name, nFolds, maxFolds, saveFile, printFlag, debugFlag, ...
    normalize);

assert(mean(PCT) > 0.5);
