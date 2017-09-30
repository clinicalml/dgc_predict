
clear all;
close all;
InitRand();

model = 'fa_lrtc';
T = GetTensor('tsize/large/large');

args = GetArgs(model, [], [], size(T));
[T_model, time, model_out] = CompleteTensor(T, model, args);

outDir = [DataDir() '/results/tsize/large/'];
save([outDir model '_final'], 'time', 'model_out', 'T_model');
