function val = J(f,u1,u2,kt)

% Evaluating the functional.
[a,  b, m, h, x, y, Nt, t, dt, nu] = parameters(kt);


for k = 1:Nt
    t1 = t(k);
    tr1 = traj1(t1);
    tr2 = traj2(t1);
    
    for i = 1:m+1
        for j = 1:m+1
            f1(i,j) = V1(x(i),y(j),tr1,tr2)*f(i,j,k) ...
                     + (nu/2)*u1(i,j,k)^2+(nu/2)*u2(i,j,k)^2;
        end
    end
    integ1(k) = sum(sum(f1)*h^2);
end

t1 = t(Nt);
tr1 = traj1(t1);
tr2 = traj2(t1);

for i = 1:m+1
    for j = 1:m+1
        f2(i,j) = V1(x(i),y(j),tr1,tr2)*f(i,j,Nt);
    end
end
integ2 = sum(sum(f2)*h^2);

val = sum(integ1)*dt + integ2;


