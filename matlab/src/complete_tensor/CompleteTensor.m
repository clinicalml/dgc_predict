function [T_model, time, model_out] = CompleteTensor(T, model, args, ...
    printFlag, debugFlag, normalize)
% model_out is any sort of output from the actual tensor completion call
% T_model is the completed tensor
% time is the runtime in seconds


if ~exist('printFlag')
    printFlag = true;
end

if ~exist('debugFlag')
    debugFlag = false;
end

if ~exist('normalize')
    normalize = false;
end

if debugFlag
    str = 'Dummy';
else
    switch(model)
        case 'tmac'
            str = 'Tmac';
        case 'asmatrix'
            str = 'AsMatrix';
        case 'constrained'
            str = 'Constrained';
        case 'mixture'
            str = 'Mixture';
        case 'ha_lrtc'
            str = 'HaLRTC';
        case 'si_lrtc'
            str = 'SiLRTC';
        case 'fa_lrtc'
            str = 'FaLRTC';
        case 'mean'
            str = 'Mean';
        case 'mean2'
            str = 'Mean2D';
        case 'matrixcomp'
            str = 'MatrixComp';
        case 'knnd'
            str = 'KNNDrug';
        case 'knnd_ts';
            str = 'KNNDrug';
        case 'knnc'
            str = 'KNNCell';
        case 'knnc_ts';
            str = 'KNNCell';
        case 'knndc'
            str = 'KNNDrugCell';
        case 'unfoldmc'
            str = 'UnfoldMC';
        case 'inftucker'
            str = 'InfTucker';
        otherwise
            error('unexpected model')
    end
end
expr = ['[T_model, model_out] = CompleteTensor' str '(T, args);'];

tic; 
if(printFlag)
    eval(expr);
else
    evalc(expr);
end
time = toc;

if normalize
    T_model = NormSigs(T_model);
end

end

