T = GetTensor('small_g978');

n = 5;
sigs = RandSigs(T, n, false);
assert(all(size(sigs) == [n,size(T,2)]));
assert(isempty(find(isnan(sigs))));

for i = 1:n
   assert(abs(norm(sigs(i,:)) - 1) < 1e-3); 
end