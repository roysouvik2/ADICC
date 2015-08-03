function [a,  b, m, h, x, y, Nt, t, dt, nu] = parameters(kt)

% Domain definition.
a = 2; 
b = 12;

% Number of subintervals
m = 60;

% Generating mesh
x = linspace(a,b,m+1);
y = linspace(a,b,m+1);

% Mesh size 
h = x(2)-x(1);

% Number of time steps in a time window.
Nt = 30;

% Final time of a time window
T = pi;

% Time grid of a time window
t = linspace(0,T,Nt);

% Time stepping.
dt = t(2)-t(1);

% Actual time grid of the simulation.
t = t+(kt-1)*T;

% Value of control parameter nu.
nu = 0.01;