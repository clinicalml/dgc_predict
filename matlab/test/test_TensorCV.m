
model = 'tmac';
tensor_name = 'small_g50';
nFolds = 132;
maxFolds = 1;
saveFile = false;
printFlag = false;
debug = false;
normalize = false;

[abs_err, cos_err] = TensorCV(model, tensor_name, nFolds, ...
    maxFolds, saveFile, printFlag, debug);

assert(abs(abs_err-.6316) < 1e-3);
assert(abs(mean(cos_err) - .4997) < 1e-3);

% tensor_name = 'allDrugs_minCellCount3';
% nFolds = 895; %6;
% debug = true;
% [abs_err1, cos_err1,~,~, testIdx1, trainIdx1, norm_true1] = ...
%     TensorCV(model, tensor_name, nFolds, maxFolds, saveFile, printFlag, debug);
% 
% tensor_name = 'allDrugs_minCellCount3_plusKD9';
% [abs_err2, cos_err2,~,~, testIdx2, trainIdx2, norm_true2] = ...
%     TensorCV(model, tensor_name, nFolds, maxFolds, saveFile, printFlag, debug);
% 
% assert(abs(mean(squeeze(norm_true1)) - 1) < 1e-3);
% assert(abs(mean(squeeze(norm_true2)) - 1) < 1e-3);
% 
% [I1,J1,K1] = ind2sub([1517 978 11], testIdx1{1});
% [I2,J2,K2] = ind2sub([2784 978 11], testIdx2{1});
% assert(isequal(I1,I2))
% assert(isequal(J1,J2))
% assert(isequal(K1,K2))
% 
% assert(isequal(abs_err1,abs_err2));
% assert(isequal(cos_err1,cos_err2));
% 
% % assert that trainIdx1 is strict subset of trainIdx2
% [I1,J1,K1] = ind2sub([1517 978 11], trainIdx1{1});
% [I2,J2,K2] = ind2sub([2784 978 11], trainIdx2{1});
% n = length(I1);
% assert(isequal(I1,I2(1:n)));
% assert(isequal(J1,J2(1:n)));
% assert(isequal(K1,K2(1:n)));




