% TEST_COMPARE - Computes and plots Figure 1
%
% See also
%  PLOT_TENSORWORKSHOP10
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt


nrep=10;
nsample = 1;
sz=[50 50 20];
trfrac=0.05:0.05:0.95;

tol=1e-3;

methods = {'matrix','constraint','mixture','tucker','tuckertrue'};
base_err = cumsum([0, 3, 1, 1, 1]);

for ll=1:nsample
  % dtrue=round(rand(1,3)*40);
  dtrue=[7,8,9];


  for kk=1:nrep
    X0=randtensor3(sz,dtrue);
    nn=prod(sz);

    for ii=1:length(trfrac)
      ntr=round(nn*trfrac(ii));
      ind=randperm(nn); ind=ind(1:ntr)';
      ind_test=setdiff(1:prod(sz), ind);
      [I,J,K]=ind2sub(sz,ind);
      yy=X0(ind);

      for mm=1:length(methods)
        switch(methods{mm})
         case 'matrix'
          %% Tensor as a matrix
          tic;
          [X1,Z1,fval1,gval1]=tensor_as_matrix(zeros(sz), {I,J,K}, yy, 1, tol);
          time(kk,ii,mm)=toc;
          gval(kk,ii,mm)=max(gval1);
          for jj=1:ndims(X0)
            Xjj = flatten_adj(Z1{jj},sz,jj);
            err(kk,ii,base_err(mm)+jj)=norm(Xjj(ind_test)-X0(ind_test));
          end
         case 'constraint'
          %% Constrained
          tic;
          [X,Z,Y,fval,gvals]=tensorconst_adm(zeros(sz),{I,J,K},yy,0,1, tol);
          time(kk,ii,mm)=toc;
          gval(kk,ii,mm)=gvals(end);
          err(kk,ii,base_err(mm)+1)=norm(X(ind_test)-X0(ind_test));
         case 'mixture'
          %% Mixture
          tic;
          [X2,Z2,fval2,gval2]=tensormix_adm(zeros(sz), {I,J,K}, yy, ...
                                            0, 1, tol);
          time(kk,ii,mm)=toc;
          gval(kk,ii,mm)=gval2(end);
          err(kk,ii,base_err(mm)+1)=norm(X2(ind_test)-X0(ind_test));
         case {'tucker','tuckertrue'}
          %% Tucker
          Xobs=zeros(sz);
          Xobs(ind)=X0(ind);
          Xobs(ind_test)=nan;
          Options(5)=100;
          if strcmp(methods{mm},'tuckertrue')
            dd = dtrue;
          else
            dd = round(dtrue*1.2);
          end
          tic;
          [Factors,G,ExplX,Xm]=tucker(Xobs, dd, Options);
          time(kk,ii,mm)=toc;
          gval(kk,ii,mm)=nan;
          err(kk,ii,base_err(mm)+1)=norm(Xm(ind_test)-X0(ind_test));
        otherwise
         error('Method [%s] unknown!', methods{mm});
        end
      end
     fprintf('frac=%g\nerr1=%s  err2=%g  err3=%g  err4=%g err5=%g\n',...
              trfrac(ii), printvec(err(kk,ii,1:3)),...
              err(kk,ii,4), err(kk,ii,5), err(kk,ii,6), err(kk,ii,7));
      fprintf('time1=%g      time2=%g time3=%g time4=%g time5=%g\n', time(kk,ii,1),time(kk,ii,2),time(kk,ii,3),time(kk,ii,4),time(kk,ii,5));
      end
  end

  file_save=sprintf('result_compare5_new_%d_%d_%d_%d_%d_%d.mat',sz(1),sz(2),sz(3),dtrue(1),dtrue(2),dtrue(3));


  save(file_save,'nrep', 'sz', 'dtrue', 'methods','err', 'trfrac','gval','time');
   
    
end

plot_tensorworkshop10
