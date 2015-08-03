function val = grad_j(u1,u2,f,p,c,kt) % Evaluating gradient of J^hat at u.

[a,b, m, h, x, y, Nt, t, dt, nu] = parameters(kt);

% Evaluating the gradient using forward difference for derivative of adjoint.

val = zeros(m+1,m+1,Nt);

for k = 1:Nt    
    for i = 2:m
        for j = 2:m
            if(c==1)
                val(i,j,k) = f(i,j,k)*(nu*(u1(i,j,k)) - (p(i+1,j,k)-p(i,j,k))/(h));
            else
                val(i,j,k) = f(i,j,k)*(nu*(u2(i,j,k)) - (p(i,j+1,k)-p(i,j,k))/(h));
            end
        end
    end
end












