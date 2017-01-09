load phege_bbs;

phe_class = 10;
geno_class = 10;


diary('phege_10_10.txt');
diary on;


g2g_meta = (g2g_seq_bbs + g2g_ppi_bbs + g2g_go_bbs + g2g_ge_bbs)/4;
omega0 = 1;
pi0 = [0.25; 0.25; 0.25; 0.25];

U0 = initialize_SNMF(rand(size(d2d_bbs, 1), phe_class),d2d_bbs,200);
V0 = initialize_SNMF(rand(size(g2g_meta, 1), geno_class),g2g_meta,200);

P{1} = d2d_bbs;
G{1, 1} = g2g_seq_bbs;
G{2, 1} = g2g_ppi_bbs;
G{3, 1} = g2g_go_bbs;
G{4, 1} = g2g_ge_bbs;
R = association;

Theta0 = rand(size(R));
Lambda0 = rand(phe_class, geno_class);

max_iter = 100;
tol = 10-4;

paramLambda = [0.01, 0.1, 1, 10];
paramdelta = [1000, 100, 10000];


t = 1;
for i = 1:4
    lambdaU = paramLambda(i);
    for j=1:4
        lambdaV = paramLambda(j);
        for p=1:3
            deltaOmega = paramdelta(p);
            for q=1:3
                deltaPi = paramdelta(q);
                
                fprintf('Experiment=%d, U = %f, V = %f, Omega = %f, Pi = %f\n', t, lambdaU,lambdaV,deltaOmega,deltaPi);
                
                tic;
                [U,V,Lambda,Theta,omega,pi, loss] = phege(U0,V0,Lambda0,omega0,pi0,Theta0,P,G,R,lambdaU,lambdaV,deltaOmega,deltaPi,max_iter,tol);
                
                result(t).runningtime =  toc;
                % settings
                result(t).phe_class = phe_class;
                result(t).geno_class = geno_class;
                result(t).lambdaU = lambdaU;
                result(t).lambdaV = lambdaV;
                result(t).deltaOmega = deltaOmega;
                result(t).deltaPi = deltaPi;
               % results
                result(t).U=U; 
                result(t).V = V;
                result(t).Lambda = Lambda;
                result(t).Theta = Theta;
                result(t).omega = omega;
                result(t).pi = pi;
                result(t).loss = loss;

                t = t +1;

            end
       end
    end
end

diary off;

save('phege_10_10.mat', 'result', '-v7.3');
