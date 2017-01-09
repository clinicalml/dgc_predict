function X = initialize_SNMF(X0,G,maxiter)
% initialize a matrix with symmetric NMF on G
X = X0;
for i=1:maxiter
    X = X.*(.5+(G*X)./(2*X*X'*X));
    %fprintf('iteration = %d, loss = %f\n',i,norm(G-X*X','fro'));
end