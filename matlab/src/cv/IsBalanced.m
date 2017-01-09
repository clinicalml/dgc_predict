function balanced = IsBalanced(trainIdx, minSigsPerDrug, sz)
balanced = true;
for i = 1:length(trainIdx)
    balanced = balanced && min(NumSigsFrom2D(trainIdx{i}, sz)) >= minSigsPerDrug;
end


end