function f1 = adj_fok_pl_adi2(u_1,u_2,kt)


[a, b,  m,  h,  x,  y,  N,  t,  dt] = parameters(kt);

clear u1 u2 u3;
nu = 0.01;
% Terminal condition for the adjoint.

f1(:,:,:) = zeros(m+1,m+1,N);

tr1 = traj(t(N));
%tr2 = traj2(t(N));
for i = 1:m+1
    for j = 1:m+1
        f1(i,j,N) = -V1(x(i),y(j),tr1);
    end
end



% Using Euler for the first time step
k = N-1;

% Interpolating u along x-axis
for j = 1:m+1
    for i = 1:m
        u_11(i,j) = inter(u_1(:,:,k),i,j,1);
        u_21(i,j) = inter(u_2(:,:,k),i,j,1);
    end
end

% Interpolating u along y-axis
for i = 1:m+1
    for j = 1:m
        u_12(i,j) = inter(u_1(:,:,k),i,j,2);
        u_22(i,j) = inter(u_2(:,:,k),i,j,2);
    end
end 


fm1 = mat_to_vec(f1(2:m,2:m,N));
length(fm1);

d = (m-1)^2;
A = sparse(d,d);
%b = zeros(d,1);

t1  = t(k);
tr1 = traj(t1);
%tr2 = traj2(t1);
for i = 2:(m)
    for j = 2:(m)
        K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h);
        K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h) - delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j));
        K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h);
        K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h) - delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j));
        
        L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h);
        L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h) - delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j));
        L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h);
        L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h) - delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1));
        
        l = index_map(i,j,m);
        
        % Building the matrix
        A(l,l) = 1 + (dt/h)*(K_3+K_2+L_2+L_3);
        
        if(i == 2)
            A(l,l+m-1) = -(dt/h)*K_2;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_2;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_3;
            else
                A(l,l-1) = -(dt/h)*L_3;
                A(l,l+1) = -(dt/h)*L_2;
            end
        elseif(i==m)
            A(l,l-(m-1)) = -(dt/h)*K_3;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_2;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_3;
            else
                A(l,l-1) = -(dt/h)*L_3;
                A(l,l+1) = -(dt/h)*L_2;
            end
        else
            A(l,l+m-1) = -(dt/h)*K_2;
            A(l,l-(m-1)) = -(dt/h)*K_3;
            A(l,l-1) = -(dt/h)*L_3;
            A(l,l+1) = -(dt/h)*L_2;
        end
        
        % Building the right hand side
        b1(l) = fm1(l) - dt*V1(x(i),y(j),tr1) - dt*nu/2*(u_1(i,j,k)^2+u_2(i,j,k)^2); 
    end
    
end
fmn = A\b1';
length(fmn);
f1(:,:,k) = vec_to_mat(fmn);
%Boundary condition.
for i=1:m+1,
    f1(i,1,k) = 0;
    f1(i,m+1,k) = 0;
    f1(1,i,k) = 0;
    f1(m+1,i,k) = 0;
end

% figure(1)
% surf(x,y,f1(:,:,k)')
% colorbar
% %view(0,90)
% pause(1)
% % 
clear b;
% Using BDF2 for the time discretization.
for k = N-2:-1:1
    
    t1 = t(k);
    tr1 = traj1(t1);
    tr2 = traj2(t1);

    
    % Interpolating u along x-axis
    for j = 1:m+1
        for i = 1:m
            u_11(i,j) = inter(u_1(:,:,k),i,j,1);
            u_21(i,j) = inter(u_2(:,:,k),i,j,1);
        end
    end
    
    % Interpolating u along y-axis
    for i = 1:m+1
        for j = 1:m
            u_12(i,j) = inter(u_1(:,:,k),i,j,2);
            u_22(i,j) = inter(u_2(:,:,k),i,j,2);
        end
    end
    
    fm1 = mat_to_vec(f1(2:m,2:m,k+1));
    fm2 = mat_to_vec(f1(2:m,2:m,k+2));
    d = (m-1)^2;
    
    A = sparse(d,d);
    %b = zeros(d,1);
    for i = 2:m
        for j = 2:m
             K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h);
             K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h) - delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j));
             K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h);
             K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h) - delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j));
             
             L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h);
             L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h) - delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j));
             L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h);
             L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h) - delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1));
             
            l = index_map(i,j,m);
            
            % Building the matrix 
            A(l,l) = 3 + 2*(dt/h)*(K_3+K_2+L_2+L_3);
            
            if(i == 2)
                A(l,l+m-1) = -2*(dt/h)*K_2;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_2;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_3;
                else
                    A(l,l-1) = -2*(dt/h)*L_3;
                    A(l,l+1) = -2*(dt/h)*L_2;
                end   
            elseif(i==m)
                A(l,l-(m-1)) = -2*(dt/h)*K_3;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_2;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_3;
                else
                    A(l,l-1) = -2*(dt/h)*L_3;
                    A(l,l+1) = -2*(dt/h)*L_2;
                end
            else
                 A(l,l+m-1) = -2*(dt/h)*K_2;
                 A(l,l-(m-1)) = -2*(dt/h)*K_3;
                 A(l,l-1) = -2*(dt/h)*L_3;
                 A(l,l+1) = -2*(dt/h)*L_2;                    
            end
            
            % Building the right hand side
            b1(l) = 4*fm1(l) - fm2(l) - 2*dt*V1(x(i),y(j),tr1)- 2*dt*nu/2*(u_1(i,j,k)^2+u_2(i,j,k)^2); 
        end
        
    end
    size(b1);    
    fmn = A\b1';
    f1(:,:,k) = vec_to_mat(fmn);
     %Boundary condition.
    for i=1:m+1,
        f1(i,1,k) = 0;
        f1(i,m+1,k) = 0;
        f1(1,i,k) = 0;
        f1(m+1,i,k) = 0;
    end
    
    
%     figure(1)
%     surf(x,y,f1(:,:,k)')
%     colorbar
%     %view(0,90)
%     pause(1)
end
