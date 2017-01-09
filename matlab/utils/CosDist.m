function [d,ny, nx] = CosDist(x,y)
d = 1 - dot(x,y)/(norm(x)*norm(y));
ny = norm(y);
nx = norm(x);
%fprintf('cos dist = %f, norm of true vec: %f, norm of est vec: %f\n', ...
%    d, norm(x), norm(y));
end
