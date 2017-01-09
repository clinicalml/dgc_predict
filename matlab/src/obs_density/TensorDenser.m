function C = TensorDenser(T, args) 

args = set_defaults(args);
C = NaN(args.numReps, 10);

for i = 1:args.numReps
    if(args.printFlag)
        disp(sprintf('starting rep #%d', i));
    end
    [test_data, test_idx, T_sub, d_vec] = SplitTensorDenser(T);
    for d = 1:length(d_vec)
            if(args.printFlag)
                disp(sprintf('  completing tensor with density %d', d_vec(d)));
            end
            model_args = GetArgs(args.model, [], [], size(T));
            T_model = CompleteTensor(T_sub{d}, args.model, model_args, ...
             args.printFlag, args.debugFlag, args.normalize);
            est = T_model(test_idx);
            C(i, d) = corr(est, test_data);
    end
end


end

function args = set_defaults(args)

    if ~isfield(args, 'numReps')
        args.numReps = 10;
    end

    if ~isfield(args, 'printFlag')
        args.printFlag = false;
    end

    if ~isfield(args, 'debugFlag')
        args.debugFlag = false;
    end

    if ~isfield(args, 'normalize')
        args.normalize = true;
    end

    if ~isfield(args, 'model')
        args.model = 'fa_lrtc';
    end
end

