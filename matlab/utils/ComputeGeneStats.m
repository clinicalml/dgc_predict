function [m,v] = ComputeGeneStats(T)

for g = 1:size(T,2)
    A = T(:,g,:);
    m(g) = nanmean(A(:));
    v(g) = nanvar(A(:));
end

end