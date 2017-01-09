InitRand();

T = GetTensor('small_g50');
args.model = 'mean2';
args.numReps = 5;

% test that each row is monotonically increasing
C = TensorDenser(T, args); 

for i = 1:size(C, 1)
    assert(issorted(C(i,1:5)));
end

% check that each column is either all NaN or all numeric
for j = 1:size(C, 2)
   assert(all(isnan(C(:,j))) || isnumeric(C(:,j))); 
end