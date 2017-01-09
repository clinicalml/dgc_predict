% PLOTRESULTS - Plots Figure 2
%
% Example
%  S=ls('result_50_50_20_*.mat','-1');
%  files=split(S,char(10));
%  [rs, frac, dim]=plotresults(files);
%
% Reference
% "On the extension of trace norm to tensors"
% Ryota Tomioka, Kohei Hayashi, and Hisashi Kashima
% arXiv:1010.0789
% http://arxiv.org/abs/1010.0789
% 
% Copyright(c) 2010 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [rs, frac, dim]=plotresults(files)

ns=length(files);

rs=zeros(ns,3);
frac=zeros(ns,1);

for ii=1:ns
  S=load(files{ii});
  sz=S.sz;
  rs(ii,:)=S.dtrue;
  ix=min(find(mean(S.err)<0.01));
  frac(ii)=S.trfrac(ix);
end
% $$$ 
% $$$ figure;
% $$$ plot(sum(rs,2),frac,'x','linewidth',2);
% $$$ for ii=1:ns
% $$$   text(sum(rs(ii,:)),frac(ii),...
% $$$        sprintf('[%d %d %d]',rs(ii,1),rs(ii,2),rs(ii,3)));
% $$$ end


rst = [min(rs(:,1), rs(:,2).*rs(:,3)),...
       min(rs(:,2), rs(:,1).*rs(:,3)),...
       min(rs(:,3), rs(:,1).*rs(:,2))];

dim=sum(rst,2);
% figure;
plot(dim,frac,'x','linewidth',2);
for ii=1:ns
  text(dim(ii),frac(ii),...
       sprintf('[%d %d %d]',rs(ii,1),rs(ii,2),rs(ii,3)));
end
p=polyfit(dim,frac,1)
hold on;
plot(xlim,polyval(p,xlim),'--','color',[.5 .5 .5],'linewidth', 2)
h=get(gca,'children');
set(gca,'children',h([2:end,1]));




%dim = prod(rs,2)+sum(rs.*(ones(ns,1)*sz),2)-sum(rs.^2,2);
%x = dim.^(1/3);
% $$$ dim(:,1) = max([rst(:,1)*(sz(1)+prod(sz(2:3)))   - rst(:,1).^2,...
% $$$            rst(:,2)*(sz(2)+prod(sz([1,3]))) - rst(:,2).^2,...
% $$$           rst(:,3)*(sz(3)+prod(sz(1:2)))   - rst(:,3).^2],[],2);
% $$$ 
% $$$ dim(:,2) = rst(:,1)*sz(1)+rst(:,2)*sz(2)+rst(:,3)*sz(3); +rst(:,1).*rst(:,2).*rst(:,3)-rst(:,1).^2-rst(:,2).^2-rst(:,3).^2;
% $$$ 
% $$$ 
% $$$ dim=min(dim,[],2);
% $$$ 
% $$$ figure, plot(dim, frac,'x', 'linewidth',2)
% $$$ for ii=1:ns
% $$$   text(dim(ii),frac(ii),...
% $$$        sprintf('[%d %d %d]',rs(ii,1),rs(ii,2),rs(ii,3)));
% $$$ end
