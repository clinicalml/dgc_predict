
rank = 3;
T = RandTensor([20 20 20], rank);

Options = [1e-4, 1, 0, 2, NaN, 1000];
err = zeros(rank,1);

for nFactors = 1:rank
    % Using evalc to suppress output.
    expression = '[~,~,err(nFactors),corcon(nFactors)]=parafac(T, nFactors, Options);';
    evalc(expression);
end

assert(err(end) < 1e-3);
assert(all(err(1:end-1) > 1));
assert(all(corcon > 99));

