
n = 10;

InitRand();
a = randperm(n);
A = rand(n);

b = randperm(n);
B = rand(n);

assert(~isequal(a,b));
assert(~isequal(A,B));

InitRand();
c = randperm(n);
C = rand(n);

assert(isequal(a,c));
assert(isequal(A,C));