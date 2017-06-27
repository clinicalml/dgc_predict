function [U, M] = tentp_probit_vem_RH(data, model, optimizer, uparam)

%% variable initialization
nvec = size(data.Y);
n = prod(nvec);
d = model.dim;
nmod = length(nvec);

% initialization
U = model.U;
M =  model.M;

iter_cut = ceil(model.em_total / 10);

for iter = 1:model.em_total
    
    if model.verbose && mod(iter,iter_cut) == 0
        fprintf('%g%% done, %d/%d iterations\n', ...
            floor(iter/model.em_total*100), iter, model.em_total);
    end
    
    % E-step
    kernel = ten_ker_func(U, uparam);
    M = optX_tentp_var(M, data.Y,  kernel, data.beta, model, uparam);
    
    tmp_Sigma0p5 = cell(nmod,1);
    inv_eigens = cell(nmod,1);
    for k = 1: nmod
        tmp_Sigma0p5{k} = kernel.V{k}* diag(double(kernel.D{k}).^(-0.5))*(kernel.V{k})';
        inv_eigens{k} = 1./kernel.D{k};
    end
    
    % M-step
    for k = 1:nmod 
        ind_other = setdiff(1:nmod, k);
        Utimsnok = ttm(M,tmp_Sigma0p5(ind_other'), ind_other');
        funObj = @(U_k)ten_gradU(U_k, k, Utimsnok, kernel.D,  kernel.V{k}, uparam, model.process);
        
        U_knew = optimizer.algorithm(funObj,U{k}(:), model.lambda* ones(nvec(k)*d,1), optimizer.options);
        
        U{k} = reshape(U_knew, nvec(k), d);
        kernel_k = ker_func(U{k}, uparam);
        kernel.V{k} = kernel_k.V;
        kernel.D{k} = kernel_k.D;
        tmp_Sigma0p5{k} = kernel.V{k}* diag(kernel.D{k}.^(-0.5))*(kernel.V{k})';
    end
    
end

% final AUC
kernel = ten_ker_func(U, uparam);
M = optX_tentp_var(M, data.Y, kernel, data.beta, model, uparam);

end