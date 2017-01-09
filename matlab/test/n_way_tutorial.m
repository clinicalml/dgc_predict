load claus

% subplot(2,1,1), h1 = plot(X(:,:,1)');
% subplot(2,1,2), h2 = plot([X(:,:,1) X(:,:,2)]');
% 
% legend(h1)
% legend(h2)
% 
% subplot(2,1,1),plot(squeeze(X(:,1,:))') 
% subplot(2,1,2),plot(squeeze(X(:,30,:))') 
% 
% clf
% subplot(3,2,1),plot(squeeze(X(2,:,:))) 
% subplot(3,2,2),plot(reshape(X(2,:),201,61)) 
% subplot(3,2,3),plot(reshape(X(2,:),201,61)') 
% subplot(3,2,4),mesh(reshape(X(2,:),201,61)) 
% subplot(3,2,5),mesh(reshape(X(2,:),201,61)')
% subplot(3,2,6),mesh(EmAx,ExAx,reshape(X(2,:),201,61)')
% axis tight


Options = [1e-6,0,1,0,10,500];

%for i = 1:10
%   [Factors{i}, it(i), err(i), cor(i)] = parafac(X,i,Options);
%end

 [ssX,Corco,It] = pftest(3,X,5,Options);
