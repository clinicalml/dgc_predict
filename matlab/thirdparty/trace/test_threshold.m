% TEST_THRESHOLD - Computes and plots Figure 2
%
% See also
%  PLOTRESULTS
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
nsample = 10;
sz=[50 50 20];
trfrac=0.05:0.05:0.95;


for ll=1:nsample
  dtrue=round(rand(1,3).*sz);
  
  err=zeros(nrep, length(trfrac));
  gval=zeros(nrep, length(trfrac));

  for kk=1:nrep

    X0=randtensor3(sz,dtrue);
    nn=prod(size(X0));


    for ii=1:length(trfrac)
      ntr=round(nn*trfrac(ii));
      ind=randperm(nn); ind=ind(1:ntr)';
      [I,J,K]=ind2sub(size(X0),ind);
      
      [X,Z,Y,fval,gvals]=tensorconst_adm(zeros(size(X0)),{I,J,K},X0(ind),0,1);

      err(kk,ii)=norm(X(:)-X0(:));
      gval(kk,ii)=gvals(end);
      fprintf('frac=%g norm(X-X0)=%g\n', trfrac(ii), err(kk,ii));
    end
  end

  file_save=sprintf('result_%d_%d_%d_%d_%d_%d.mat',sz(1),sz(2),sz(3),dtrue(1),dtrue(2),dtrue(3));

  save(file_save,'nrep', 'sz', 'dtrue', 'err', 'trfrac','gval');

end

S=ls('result_50_50_20_*.mat','-1');
files=split(S,char(10));
[rs, frac, dim]=plotresults(files);
