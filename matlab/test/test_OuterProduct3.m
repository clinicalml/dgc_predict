x = 1;
y = [1 2];
z = [1 2 3];

T = OuterProduct3(x,y,z);

assert(isequal(squeeze(T(1,1,:)), [1 2 3]'));
assert(isequal(squeeze(T(1,2,:)), [2 4 6]'));