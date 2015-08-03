% model definition

% griglia di x1 e x2. D1 e D2 sono i punti dei lati del dominio
[m.xx1, m.xx2, m.tt] = meshgrid(D2,D1,TT);

% first component
m.b1 = @(x1,x2,u1,u2,t) interp3(m.xx1, m.xx2, m.tt, u1, x2, x1,t)  ;  % interpola la superficie del controllo

m.sig11 = @(x1,x2)  1;
m.sig12 = @(x1,x2)  0;

m.dsig11 = @(x1,x2)  0;
m.dsig12 = @(x1,x2)  0;


%second component
m.b2 = @(x1,x2,u1,u2,t) interp3(m.xx1, m.xx2, m.tt, u2, x2, x1,t)  ;  % interpola la superficie del controllo

m.sig21 = @(x1,x2)  0;
m.sig22 = @(x1,x2)  1;

m.dsig21 = @(x1,x2)  0;
m.dsig22 = @(x1,x2)  0;