function f1 = fok_pl(u_1,u_2,fin,kt)

% Loading parameters.
[a, b,  m,  h,  x,  y,  N,  t,  dt, nu] = parameters(kt);

clear u1 u2 u3;

% Initial condition.
f1(:,:,1) = fin;

% Using Euler for the first time step
k = 2;

% Mapping the 2D solution matrix to 1D solution vector.
fm1 = mat_to_vec(f1(2:m,2:m,1));

% Number of grid points along each side is m+1.
% Number of boundary points along each side is 2.
% Hence total number of degrees of freedom along each side is m-1.
d = (m-1)^2;

A = sparse(d,d);

t1  = t(k);

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
        
        % Map of matrix index to vector index.
        l = index_map(i,j,m);
        
        % Building the matrix
        A(l,l) = 1 + (dt/h)*(K_3+K_2+L_2+L_3);
        
        if(i == 2)
            A(l,l+m-1) = -(dt/h)*K_1;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_1;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_4;
            else
                A(l,l-1) = -(dt/h)*L_4;
                A(l,l+1) = -(dt/h)*L_1;
            end
        elseif(i==m)
            A(l,l-(m-1)) = -(dt/h)*K_4;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_1;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_4;
            else
                A(l,l-1) = -(dt/h)*L_4;
                A(l,l+1) = -(dt/h)*L_1;
            end
        else
            A(l,l+m-1) = -(dt/h)*K_1;
            A(l,l-(m-1)) = -(dt/h)*K_4;
            A(l,l-1) = -(dt/h)*L_4;
            A(l,l+1) = -(dt/h)*L_1;
        end
        
        % Building the right hand side
        b1(l) = fm1(l);
    end
    
end

% Calculating the solution.
fmn = A\b1';

% Mapping the 1D solution vector to the original 2D solution matrix
f1(:,:,k) = vec_to_mat(fmn);

%Boundary condition.
for i=1:m+1,
    f1(i,1,k) = 0;
    f1(i,m+1,k) = 0;
    f1(1,i,k) = 0;
    f1(m+1,i,k) = 0;
end

% Plotting solutions and evaluating integral.

%sum(trapz(f1(:,:,k)))*h^2
% figure(1)
% surf(x,y,f1(:,:,k)')
% colorbar
% view(0,90)
% pause(1)


% Using BDF2 for the time discretization for subsequent steps
for k = 3:N
    
    t1 = t(k);
    
    
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
    
    % Mapping the 2D solution matrix to 1D solution vector.
    fm1 = mat_to_vec(f1(2:m,2:m,k-1));
    fm2 = mat_to_vec(f1(2:m,2:m,k-2));
    d = (m-1)^2;
    
    A = sparse(d,d);
    
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
                A(l,l+m-1) = -2*(dt/h)*K_1;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_1;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_4;
                else
                    A(l,l-1) = -2*(dt/h)*L_4;
                    A(l,l+1) = -2*(dt/h)*L_1;
                end
            elseif(i==m)
                A(l,l-(m-1)) = -2*(dt/h)*K_4;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_1;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_4;
                else
                    A(l,l-1) = -2*(dt/h)*L_4;
                    A(l,l+1) = -2*(dt/h)*L_1;
                end
            else
                A(l,l+m-1) = -2*(dt/h)*K_1;
                A(l,l-(m-1)) = -2*(dt/h)*K_4;
                A(l,l-1) = -2*(dt/h)*L_4;
                A(l,l+1) = -2*(dt/h)*L_1;
            end
            
            % Building the right hand side
            b1(l) = 4*fm1(l) - fm2(l);
        end
        
    end
    
    % Calculating the solution.
    fmn = A\b1';
    
    % Mapping the 1D solution vector to the original 2D solution matrix.
    f1(:,:,k) = vec_to_mat(fmn);
    
    %Boundary condition.
    for i=1:m+1,
        f1(i,1,k) = 0;
        f1(i,m+1,k) = 0;
        f1(1,i,k) = 0;
        f1(m+1,i,k) = 0;
    end
    
    % Plotting solutions and evaluating integral.
    
    %     abs(1-sum(trapz(f1(:,:,k)))*h^2)
    %     figure(2)
    %     surf(x,y,f1(:,:,k)')
    %     colorbar
    %     view(0,90)
    %     pause(1)
end
