%model = {'tmac', 'asmatrix', 'constrained', 'ha_lrtc', 'fa_lrtc', 'si_lrtc'};
nFolds = 10;
maxFolds = 10;
saveFile = false;
printFlag = true;
debugFlag = false;
normalize = true;

T = GetTensor('tsize/small/small');

%$T = T(1:10,1:20,1:3);

[T_imputed, PCT, PCTf, time, testIdx] = TensorCV(model, T, 'benchmarking', nFolds, ...
                                          maxFolds, saveFile, printFlag, debugFlag, normalize)

save(DataDir(sprintf('results/tsize/small/benchmarking/%s_full_small_tensor.mat', model)), 'model', 'PCT','PCTf','time','testIdx');
