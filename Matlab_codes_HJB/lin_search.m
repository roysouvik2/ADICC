function alpha = lin_search(f,u1,u2,d1,d2,g1,g2,rho,kt)

% Performing the linear search algorithm.
[a,b, m, h, x, y, N, t, dt,nu] = parameters(kt);

% Starting value of alpha.
alpha = 0.5;

delta = 0.1 ;

s = 1;
s_max = 10;

while( s < s_max)
    fn = fok_pl(u1+alpha*d1, u2+alpha*d2,rho,kt);
    
    c1 = J(fn,u1 + alpha*d1,u2+alpha*d2,kt);
    c2 = J(f,u1,u2,kt) + delta*alpha*inner_g(g1,g2,d1,d2);
    
    % Checking for Armijo criterion.
    if(c1 <= c2 )
        return; 
    else
        alpha=alpha/2;
    end
    s = s + 1;
end
if(s == s_max)
    alpha = 0;
end
