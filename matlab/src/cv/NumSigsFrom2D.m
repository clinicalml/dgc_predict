function [numSigs] = NumSigsFrom2D(idx2D, tensorSize)

[I, ~] = ind2sub([tensorSize(1) tensorSize(3)], idx2D);

numSigs = zeros(tensorSize(1), 1);

for drug = 1:tensorSize(1)
    numSigs(drug) = length(find(I == drug));
end

end