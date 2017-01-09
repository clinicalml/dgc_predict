function sim_mod = bbs_normalize(sim, Es, max_iter)
% BBS_NORMALIZE - normalize the similarities by bi-stochastic algorithm.
%
%Input:
% sim - original similarity matrix
% Es - convergence parameter. e.g., Es = 1e-6
% max_iter - maximum iteration parameter. e.g., max_iter = 200

g0 = sim;
n=size(sim,1);
iter = 0;
diff = 1;

while diff >= Es
    if iter <= max_iter
        iter = iter + 1;

        temp1=(ones(n,1)*ones(1,n)*g0)/n;
        temp2 = eye(n)-g0;

        g1 = g0+((temp1+temp2)*ones(n,1)*ones(1,n))/n-temp1;
        g1(find(g1<0)) = 0;
    else
        break
    end
    
    J1 = trace(g1'*g1 - 2*sim'*g1);
    J0 = trace(g0'*g0 - 2*sim'*g0);
    
    diff = abs(J1-J0);

    g0 = g1;
    
    disp(sprintf('iteration = %d, loss = %f ', iter, abs(J1)  ));

end

sim_mod = g0;