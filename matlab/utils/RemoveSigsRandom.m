 function [T_remaining, sigs, pairs, T_removed] = RemoveSigsRandom(...
     T, n, keep1, file, printFlag)
 % the keep1 (boolean) argument, if set to true, removes a random set of
 % signatures but ensures that there remains at least one signature per
 % drug (1st mode of tensor)
 
 T_remaining = T;
 T_removed = nan(size(T));
 
 if ~exist('keep1')
     keep1 = false;
 end
 
 if ~exist('printFlag')
     printFlag = false;
 end
 
 if ~keep1
     [~, I, J] = RandSigs(T, n);
     
     sigs = NaN(n, size(T,2));
     
     for i = 1:n
         x = T(I(i),:,J(i));
         sigs(i,:) = x;
         T_removed(I(i), :, J(i)) = x;
         T_remaining(I(i), :, J(i)) = NaN;      
     end
     
     pairs = [I J];
 else
     % check that we are starting with enough signatures to be able to keep
     % at least one per drug
     assert(NumSigs(T) - n >= size(T,1));
     
     % randomly order all signatures in tensor
     [~, I, J] = RandSigs(T, NumSigs(T));
     
     % initialize variables
     sigs = NaN(n, size(T,2));
     num_removed = 0;
     idx = 1;
     A = squeeze(~isnan(T(:,1,:)));
     drug_counts = sum(A, 2);
     pairs = NaN(n, 2);
     
     while(num_removed < n)
        % check if we can remove sig idx
        if(drug_counts(I(idx)) > 1)
            if(printFlag)
                fprintf('%d works\n', idx);
            end
            num_removed = num_removed + 1;
            drug_counts(I(idx)) = drug_counts(I(idx)) - 1;
            
            x = T(I(idx),:,J(idx));
            sigs(num_removed, :) = x;
            T_removed(I(idx), :, J(idx)) = x;
            T_remaining(I(idx),:,J(idx)) = NaN;
            
            pairs(num_removed,:) = [I(idx), J(idx)];
        else
           if(printFlag)
               fprintf('%d doesnt work\n', idx); 
           end
        end
        idx = idx + 1;
     end
 end
 

 if(exist('file'))
    base = '../data/results/allDrugs_minCellCount6/cv_loo/';
    file_sigs = [base file];
    file_pairs = [base '/drug_cell_pairs/heldout_pairs.txt'];
    dlmwrite(file_sigs, sigs, ' ');
    dlmwrite(file_pairs, pairs, ' ');
 end
 
 
 end