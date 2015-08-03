% test controlled MC 2 dim
clear all
close all

fname='data_traj_cir.mat';
load(fname);   % all data needed for the simulation



% number of desired trajectories
Nsample=15;
Ntime = 300; % number of time points for the Milstein integration


T=T1-T0;

[ Nx1, Nx2, Nt] = size(u1)

D1 = linspace(a,b,Nx1);
D2 = linspace(a,b,Nx2);
TT = linspace(0,T,Nt);

% load model definitions. Modify inside if you want
% change the stochastic model
model

x00 = ip(1);  x01=ip(2);


% Integrate the stochastic process
[y, t] = milstein_2D(m,x00,x01,u1,u2,Nsample,T,Ntime);

 plot(y(:,:,1),y(:,:,2),'Color',[0.5 0.5 0.5]);
 hold on;
 axis([a b a  b]);
 