
A = load('/Users/rhodos/Desktop/Research/LINCS/data/expr/tensor/from_jingshu/T_train_d600V2.mat');
T_train = A.out;

A = load('/Users/rhodos/Desktop/Research/LINCS/data/expr/tensor/from_jingshu/T_train_d600V2.mat');
T_test = A.out;

Tsm_train = T_train(1:10,1:10,:);
Tsm_test = T_test(1:10,1:10,:);

args = GetArgs('inftucker');

tic;
[T_est, out] = CompleteTensorInfTucker(Tsm_train, args);
toc;

