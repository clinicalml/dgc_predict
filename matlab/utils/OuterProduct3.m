function T = OuterProduct3(x,y,z)

T = nan(length(x), length(y), length(z));

for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            T(i,j,k) = x(i)*y(j)*z(k);
        end
    end
end

end