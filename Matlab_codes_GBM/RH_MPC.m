% Main function for running the MPC-algorithm.

[a,  b, m, h, x, y, Nt, t, dt, nu] = parameters(1);

% Total number of grid points along each direction.
Nx = m+1;

% Starting point of the trajectory.
ip = traj1(t(1));

% Assigning initial pdf.
rho = zeros(Nx,Nx);

for i = 2:Nx-1
    for j = 2:Nx-1
        rho(i,j) = exp((-(x(i)-ip(1))^2-(y(j)-ip(2))^2)/0.5);
    end
end
for i=1:Nx,
    rho(i,1) = 0;
    rho(i,Nx) = 0;
    rho(1,i) = 0;
    rho(Nx,i) = 0;
end

% Normalizing the integral of the initial PDF to 1.
rho = rho/(sum(sum(rho))*h^2);

% Number of time windows.
Nw = 1;

% Assigning initial guess for u.
u1 = zeros(Nx,Nx,Nw*Nt);
u2 = zeros(Nx,Nx,Nw*Nt);

% Final time in a time window.
T = t(Nt);

% Initial time.
T0 = t(1);
% Final time of simulation.
T1 = t(Nt)+(Nw-1)*T;

% Initializing PPDF of solutions.
sol = zeros(Nx,Nx,Nt); 

for step = 1:Nw
    step
    [u1_t,u2_t] = ncg(rho,step);
    
    % Appending the value of u in the current time window to the value of u
    % in the previous time window.
    u1(:,:,1+(step-1)*Nt:Nt*step) = u1_t;
    u2(:,:,1+(step-1)*Nt:Nt*step) = u2_t;
    
    sol = fok_pl(u1_t,u2_t,rho,step);
        
    % Assigning the solution at the final time step of the time window as
    % the initial condition of the next time window.
    rho = sol(:,:,Nt);
    
    
    
end
fname = 'data_traj_cir.mat';
save(fname,'a','b','Nx','T0','T1','Nt','u1','u2','ip');

    




    