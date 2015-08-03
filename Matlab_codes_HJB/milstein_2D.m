function   [y, t] = milstein_2D(M,x1_ini,x2_ini,u1,u2,Nsample,T,Ntime)
%
%   milstein_2D(M,x1_ini,x2_ini,u1,u2,Nsample,T,Ntime)
%
%   M = struct with the model definition
%   x1_ini, x2_ini = initial data 
%   u1(x1,x2,t), u2(x1,x2,t) = controls 
%   Nsample = number of samples
%   T   = final time
%   Ntime = number of time point for the Milstein integration. Ex: Ntime=100
%


rng('shuffle');   % changes the seed of the random generator

t=linspace(0,T,Ntime);
dt=t(2)-t(1);
dt2=sqrt(dt);

y = zeros(Ntime,Nsample,2);
m=1;


while m <= Nsample    % repeat for all samples
    
    
    n=1;
    
    x10=x1_ini;
    x20=x2_ini;
    
    y(n,m,1)=x10;
    y(n,m,2)=x20;  % save results 
    
     
    dW = dt2 * randn(Ntime,2); % prepare the random values
    
     % integrates with the Milstein method
    for tt=t
    
        dw = dW(n,:);    
    
       x1 = x10 + M.b1(x10,x20,u1,u2,tt) * dt + ...
                    M.sig11(x10,x20) * dw(1) + ...
                    M.sig12(x10,x20) * dw(2) + ...
                    0.5*M.sig11(x10,x20) * M.dsig11(x10,x20) * (dw(1)*dw(1)-dt) + ...
                    0.5*M.sig12(x10,x20) * M.dsig12(x10,x20) * (dw(2)*dw(2)-dt);
        
        x2 = x20 + M.b2(x10,x20,u1,u2,tt) * dt + ...
                    M.sig21(x10,x20) * dw(1) + ...
                    M.sig22(x10,x20) * dw(2) + ...
                    0.5*M.sig21(x10,x20) * M.dsig21(x10,x20) * (dw(1)*dw(1)-dt) + ...
                    0.5*M.sig22(x10,x20) * M.dsig22(x10,x20) * (dw(2)*dw(2)-dt);
    
        
        % update values   
        x10=x1;
        x20=x2;
        
        
        y(n,m,1)=x10;  % save the results     
        y(n,m,2)=x20; 
        n=n+1;
        
    end  % end Milstein
    m=m+1;

end


