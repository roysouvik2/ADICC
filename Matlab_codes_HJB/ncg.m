function [u1, u2] = ncg(rho,kt)

[a,  b, m, h, x, y, Nt, t, dt, nu] = parameters(kt);

u1 = zeros(m+1,m+1,Nt);
u2 = zeros(m+1,m+1,Nt);

f = fok_pl(u1,u2,rho,kt);
p = adj_fok_pl(u1,u2,kt);

% Evaluating the gradient with the initial guess.
g1 = grad_j(u1,u2,f,p,1,kt);
g2 = grad_j(u1,u2,f,p,2,kt);

d1 = -g1;
d2 = -g2;

fprintf('Initial norm of gradient = %.3e\n\n', sqrt(inner_g(g1,g2,g1,g2)));
fprintf('Starting the NCG\n\n');

k = 1;
k_max = 50;
tol = 10^-4;

while(k < k_max && sqrt(inner_g(g1,g2,g1,g2)) > tol)
    alpha = lin_search(f,u1,u2,d1,d2,g1,g2,rho,kt);
    if(alpha < eps)
        fprintf('Line Search terminates\n');
        fprintf('k=%d, alpha= %.4e,grad(J)=%.3e, J(u) = %.2e\n',k,alpha,sqrt(inner_g(g1,g2,g1,g2)), J(f,u1,u2,kt));
        break;
    else
        %fprintf('Line search gives alpha = %.4e\n',alpha);
        u1 = u1 + alpha*d1;
        u2 = u2 + alpha*d2;
        %inner_g(u1,u2,u1,u2)
        f = fok_pl(u1,u2,rho,kt);
        p = adj_fok_pl(u1,u2,kt);
        
        gn1 = grad_j(u1,u2,f,p,1,kt);
        gn2 = grad_j(u1,u2,f,p,2,kt);
        
        % Evaluating beta by Dai-Yuan method.
        beta = inner_g(gn1,gn2,gn1,gn2)/inner_g(d1,d2,gn1-g1,gn2-g2);
        
        % Evaluating beta by Hager-Zhang method
        %nor = inner_g(d1,d2,d1,d2)/inner_g(d1,d2,gn1-g1,gn2-g2);
        %term1 = gn1-g1 - 2*d1*nor;
        %term2 = gn2-g2 - 2*d2*nor;
        
        %beta = inner_g(term1,term2,gn1,gn2)/inner_g(d1,d2,gn1-g1,gn2-g2);
        
        % Descent directions
        d1 = -gn1 + beta*d1;
        d2 = -gn2 + beta*d2;
        
    end
    
    % Evaluating gradient for next iteration
    g1 = gn1;
    g2 = gn2;
    fprintf('k=%d, alpha= %.4e,grad(J)=%.3e, J(u) = %.2e\n',k,alpha,sqrt(inner_g(g1,g2,g1,g2)), J(f,u1,u2,kt));
    k = k + 1;
end
