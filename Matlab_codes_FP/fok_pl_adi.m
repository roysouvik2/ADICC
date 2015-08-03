% Solving the discrete forward Fokker-Planck equation.
function f1 = fok_pl_adi(u_1,u_2,fin,kt)

[a, b,  m,  h,  x,  y,  N,  t,  dt, nu] = parameters(kt);

% load 'data_lin_near.mat';
% u_1 = u1;
% u_2 = u2;

f1(:,:,1) = fin;


% Interpolating u along x-axis
for k = 1:N
    for j = 1:m+1
        for i = 1:m
            u_11(i,j,k) = inter(u_1(:,:,k),i,j,1);
            u_21(i,j,k) = inter(u_2(:,:,k),i,j,1);
        end
    end
end

% Interpolating u along y-axis
for k = 1:N
    for i = 1:m+1
        for j = 1:m
            u_12(i,j,k) = inter(u_1(:,:,k),i,j,2);
            u_22(i,j,k) = inter(u_2(:,:,k),i,j,2);
        end
    end
end


for k = 2:N
    k;
    t1 = t(k)-dt/2;
    % Interpolating at half time step.
    u_111 = (u_11(:,:,k-1)+u_11(:,:,k))/2;
    u_131 = u_12(:,:,k);
    u_121 = u_12(:,:,k-1);
    
    u_212 = (u_21(:,:,k-1)+u_21(:,:,k))/2;
    u_232 = u_22(:,:,k);
    u_222 = u_22(:,:,k-1);
       
    % Boundary condition.
    for i=1:m+1                           
        f2(i,1) = 0;
        f2(i,m+1) = 0;
        f2(1,i) = 0;
        f2(m+1,i) = 0;
    end
    % Implicit difference in x-direction
    for j = 2:m
        A = sparse(m-1,m-1);b =zeros(m-1,1);
        for i = 2:m
            K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j)))*B_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h);
            K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h) - delta_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j))*B_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j));
            K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j)))*B_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h);
            K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h) - delta_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j))*B_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j));
            
            L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_121(i,j),u_222(i,j)))*B_2(t1,x(i),y(j)+h/2,h,u_121(i,j),u_222(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h);
            L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h) - delta_2(t1,x(i),y(j)+h/2,h,u_121(i,j),u_222(i,j))*B_2(t1,x(i),y(j)+h/2,h,u_121(i,j),u_222(i,j));
            L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_121(i,j-1),u_222(i,j-1)))*B_2(t1,x(i),y(j)-h/2,h,u_121(i,j-1),u_222(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h);
            L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h) - delta_2(t1,x(i),y(j)-h/2,h,u_121(i,j-1),u_222(i,j-1))*B_2(t1,x(i),y(j)-h/2,h,u_121(i,j-1),u_222(i,j-1));
            
            b(i-1) = 2*f1(i,j,k-1)/dt + (1/h)*(L_1*f1(i,j+1,k-1)-(L_2+L_3)*f1(i,j,k-1)+L_4*f1(i,j-1,k-1));
            if(i==2)
                A(i-1,i) = -K_1/h;
                A(i-1,i-1) = 2/dt+K_2/h+K_3/h;
            elseif(i==m)
                A(i-1,i-2) = -K_4/h;
                A(i-1,i-1) = 2/dt+K_2/h+K_3/h;
            else
                A(i-1,i-1) = 2/dt+K_2/h+K_3/h;
                A(i-1,i) = -K_1/h;
                A(i-1,i-2) = -K_4/h;              
            end
        end
        
        unew = A\b; 
        for i=1:m-1,
            f2(i+1,j) = unew(i);
        end
        
    end
    
    for i=1:m+1,                              
        f1(i,1,k) = 0;
        f1(i,m+1,k) = 0;
        f1(1,i,k) = 0;
        f1(m+1,i,k) = 0;
    end
    % Implicit difference in y-direction
    for i = 2:m
        A = sparse(m-1,m-1);b =zeros(m-1,1);
        for j = 2:m
            K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j)))*B_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h);
            K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h) - delta_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j))*B_1(t1,x(i)+h/2,y(j),h,u_111(i,j),u_212(i,j));
            K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j)))*B_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h);
            K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h) - delta_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j))*B_1(t1,x(i)-h/2,y(j),h,u_111(i-1,j),u_212(i-1,j));
            
            L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_131(i,j),u_232(i,j)))*B_2(t1,x(i),y(j)+h/2,h,u_131(i,j),u_232(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h);
            L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h) - delta_2(t1,x(i),y(j)+h/2,h,u_131(i,j),u_232(i,j))*B_2(t1,x(i),y(j)+h/2,h,u_131(i,j),u_232(i,j));
            L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_131(i,j-1),u_232(i,j-1)))*B_2(t1,x(i),y(j)-h/2,h,u_131(i,j-1),u_232(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h);
            L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h) - delta_2(t1,x(i),y(j)-h/2,h,u_131(i,j-1),u_232(i,j-1))*B_2(t1,x(i),y(j)-h/2,h,u_131(i,j-1),u_232(i,j-1));
            
            b(j-1) = 2*f2(i,j)/dt + (1/h)*(K_1*f2(i+1,j)-(K_2+K_3)*f2(i,j)+K_4*f2(i-1,j));
            if(j==2)
                A(j-1,j) = -L_1/h;
                A(j-1,j-1) = 2/dt+L_2/h+L_3/h;
            elseif(j==m)
                A(j-1,j-2) = -L_4/h;
                A(j-1,j-1) = 2/dt+L_2/h+L_3/h;
            else
                A(j-1,j-1) = 2/dt+L_2/h+L_3/h;
                A(j-1,j) = -L_1/h;
                A(j-1,j-2) = -L_4/h;              
            end
        end
        
        unew = A\b;   
        for j=1:m-1,
            f1(i,j+1,k) = unew(j);
        end
        
    end

end


 

 


