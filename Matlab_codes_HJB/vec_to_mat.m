function g = vec_to_mat(f)

a = length(f);
b = sqrt(a);
g = zeros(b+2,b+2);

for i = 2:b+1
    for j = 2:b+1
        g(i,j) = f((b)*(i-2)+j-1);
    end
end
