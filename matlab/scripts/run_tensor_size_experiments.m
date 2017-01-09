
for t = 1:length(tsizes)
  tsize = tsizes{t};
  out = load(DataDir(sprintf('expr/tensor/tsize/%s/%s.mat', tsize, tsize)));
  T = out.T;
  models = {'mean', 'mean2','knnd','fa_lrtc'};
  [T_imp, PCT, PCTf, time, testIdx] = TensorCV4(models, T);

  save(DataDir(sprintf('/results/tsize/%s/%s_tensor_results.mat', tsize, tsize)) , 'PCT', 'PCTf', 'T_imp', 'time', 'testIdx','-v7.3');
end

