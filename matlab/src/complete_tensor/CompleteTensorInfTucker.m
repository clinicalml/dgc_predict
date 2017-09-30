function [T_model, out] = CompleteTensorInfTucker(T, args)

        warning('The optimization doesnt seem to be working here..')
        
        % convert T to a tensor object
        T = tensor(T);

        % setup data struct
        data = args.data;
        data.Y = T;
        
        % setup model
        model = args.model;
       
        model.U = cell(3,1);
         for j = 1:3
             fprintf('j=%d\n',j);
             % model.U{j} = nvecs(T, j, args.dim); 
             model.U{j} = randn(size(T,j),model.dim);
         end
        model.M = T;
        
        % run tensor completion
        [out, T_model] = tentp_probit_vem_RH(data, model, args.optimizer, args.uparam);
end
