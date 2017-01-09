
%args = {};
%args.model = 'mean2';
args.printFlag = true;
args.numReps = 5;

for size = 5 % 1:5
   fprintf('size = %d\n', size);
   for rep = 1:5
        if(size == 5)
           rep_str = 1;
        else
           rep_str = rep; 
        end
        fprintf('rep = %d\n', rep);
        T = GetTensor('tsize/small/small');
        C = TensorDenser(T, args);
        outfile = sprintf('%s/results/tsize/small/obs_density/%s_C_size%d_rep%d.mat', ...
            DataDir(), args.model, size, rep);
        save(outfile, 'C');
   end
end







