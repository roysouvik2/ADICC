% Checking the solution of the forward Fokker-Planck equation.


[a, b,  m,  h,  x,  y,  N,  t,  dt, nu] = parameters(1);


% Assigning initial pdf.
for i = 1:m+1
    for j = 1:m+1
        rho(i,j) = exp((-(x(i)-10)^2-(y(j)-10)^2)/0.2);
    end
end

rho = rho/(sum(sum(rho))*h^2);
figure(2)
surf(x,y,rho(:,:)')
colorbar
view(0,90)
pause(1)

for tim = 1:1
    tim
    for k = 1:N
        for i = 1:m+1
            for j = 1:m+1
                u1(i,j,k) = 1/x(i);
                u2(i,j,k) = 1/y(j);
            end
        end
    end
    f = fok_pl(u1,u2,rho,1);
    %fok_pl_direct(u1,u2,rho);
    %p1 = adj_fok_pl_adi2(u1,u2,1);

end


 




