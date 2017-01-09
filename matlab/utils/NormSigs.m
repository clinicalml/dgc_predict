function T = NormSigs(T)
% Normalizes each gene expression signature to have unit length

for i = 1:size(T,1)
  for k = 1:size(T,3)
    if abs(norm(T(i,:,k)) - 1) > 1e-3
      T(i,:,k) = T(i,:,k) / norm(T(i,:,k));
      %fprintf('normalizing (%d,%d)\n', i,k);
    end
  end
end

end
