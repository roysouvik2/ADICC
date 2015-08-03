function val = inner_g(g1,g2,d1,d2)

% Evaluating inner product of a vector function.
[a,  b, m, h, x, y, Nt, t, dt, nu] = parameters(1);


for  k = 1:Nt
    for i = 1:m+1
        for j = 1:m+1
            f1(i,j) = g1(i,j,k)*d1(i,j,k)+g2(i,j,k)*d2(i,j,k);
        end
    end
    integ(k) = sum(sum(f1)*h^2);
end
val = sum(integ)*dt;

    
       
    