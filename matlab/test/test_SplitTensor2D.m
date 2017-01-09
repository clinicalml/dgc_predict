% This also tests ConstructTrainingTensorFrom2D and NumSigsFrom2D.

printFlag = false;

T = GetTensor('d6_subset');


% First test that SplitTensor and SplitTensor2D give the same results in a
% case where the split is balanced on the first try

InitRand();
[test1, train1] = SplitTensor(T);

InitRand();
[test2, train2] = SplitTensor2D(T);

T1 = ConstructTrainingTensor(T, train1{1});
T2 = ConstructTrainingTensorFrom2D(T, train2{1});

assert(isequaln(T1, T2));

% Test NumSigsFrom2D

numSigs1 = NumSigs(T1, 'drug');
numSigs2 = NumSigsFrom2D(train2{1}, size(T));

assert(isequal(numSigs1, numSigs2));

% Now find cases where SplitTensor2D has to reshuffle, and verify that
% SplitTensor returned a result that was not balanced, while
% SplitTensor2D's result is balanced

for minSigsPerDrug = 1:2

for i = 60:63
    InitRand(i)
    T_sparse = RemoveSigsRandom(T, floor(NumSigs(T)/3));
    
    if min(NumSigs(T_sparse,'drug')) > minSigsPerDrug
        if printFlag fprintf('i=%d\n',i); end
        
        InitRand();
        [test1, train1] = SplitTensor(T_sparse);

        InitRand();
        [test2, train2, nIter] = SplitTensor2D(T_sparse, 10, 10, minSigsPerDrug);
        
        assert(IsBalanced(train2, minSigsPerDrug, size(T)));
        
        if nIter == 1
            
            % check that the resulting tensors are equal and balanced
            for j = 1:10
                T1 = ConstructTrainingTensor(T_sparse, train1{j});
                T2 = ConstructTrainingTensorFrom2D(T_sparse, train2{j});
                assert(isequaln(T1, T2));
                assert(min(NumSigs(T2,'drug')) >= minSigsPerDrug);
            end
        else
            % check that the resulting tensors are unequal and that the
            % former is imbalanced while the latter is balanced
            m = Inf;
            for fold = 1:length(train1)
                T1 = ConstructTrainingTensor(T_sparse, train1{fold});
                T2 = ConstructTrainingTensorFrom2D(T_sparse, train2{fold});
                assert(min(NumSigs(T2,'drug')) >= minSigsPerDrug);
                m = min(m, min(NumSigs(T1,'drug')));
            end
            assert(m < minSigsPerDrug);
        end
    end
    
end
    
end





