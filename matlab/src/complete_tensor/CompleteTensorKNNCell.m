function [T, out] = CompleteTensorKNNCell(T, args)

%% Compute cell-cell correlation matrix using all data in tensor
%M = Unfold(T, size(T), 3)';
%C = corr(M, 'rows', 'pairwise');

% % This way is faster and gives almost identical results
% C = nan(size(T,3),size(T,3),size(T,1)); 
% for d = 1:size(T,1)
%     C(:,:,d) = corr(squeeze(T(d,:,:))); 
% end
% C = mean(C, 3, 'omitnan');
% C_save = C;

% This is faster and is also memory efficient
C = zeros(size(T,3), size(T,3));
nanCt = zeros(size(T,3), size(T,3));
for d = 1:size(T,1)
    newCors = corr(squeeze(T(d,:,:)));
    tmp = cat(3,C,newCors);
    C = nansum(tmp,3);
    nanCt = nanCt + isnan(newCors);
end
C = C ./ (size(T, 1) - nanCt);
clear nanCt tmp newCors

% set any negative correlations to 0 so that they don't contribute 
if ~isempty(find(C(:) <= 0, 1))
   fprintf('Setting negative correlations to 0.\n');
   C(C < 0) = 0;
end

% set diagonal to 0
C(eye(size(C))~=0) = 0;

% set NaN's to 0
C(isnan(C)) = 0;

%%% Using this correlation matrix as the similarity measure, estimate each
%%% missing profile as a weighted combination of profiles for the same drug
%%% in other cell types, using the K nearest cell types

for drug = 1:size(T,1)
    cellsMissing = find(isnan(T(drug,1,:)));
        
    K = min(args.K, size(T,3) - length(cellsMissing));
    if K < args.K
        fprintf('Warning: only %d profiles available for drug %d\n', K, drug);
    end
    
    for c = 1:length(cellsMissing)
        cell = cellsMissing(c);
        preweights = C(cell,:);
        preweights(cellsMissing) = 0;
        [preweights_sort, idx_sort] = sort(preweights, 'descend');
        weights = preweights_sort(1:K) / sum(preweights_sort(1:K));
        T(drug,:,cell) = weights * squeeze(T(drug,:,idx_sort(1:K)))';
    end
end

%%% If any missing elements still remain, use Mean-1D
if length(find(isnan(T))) > 0
    T = CompleteTensorMean(T, args);
end

out = [];

end