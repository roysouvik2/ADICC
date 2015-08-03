function f = mat_to_vec(g)

[a,b] = size(g);


count = 1;

for i = 1:a
    for j = 1:b
        f(count) = g(i,j);
        count = count+1;
    end
end
