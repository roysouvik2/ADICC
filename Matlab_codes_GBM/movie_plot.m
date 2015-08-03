function val = movie_plot(u1,u2)

[a,  b, m, h, x, y, N, t, dt] = parameters(1);
rho = zeros(m+1,m+1);
% Assigning initial pdf.
for i = 2:m
    for j = 2:m
        rho(i,j) = exp((-(x(i))^2-(y(j))^2)/0.5);
    end
end
for i=1:m+1,
    rho(i,1) = 0;
    rho(i,m+1) = 0;
    rho(1,i) = 0;
    rho(m+1,i) = 0;
end
rho = rho/(sum(sum(rho))*h^2);

% Velocity Vectors.
vid = VideoWriter('circle.avi');
vid.FrameRate = 5;
open(vid);

figure(1)
F(N-1) = struct('cdata',[],'colormap',[]);

for k = 1:N
    quiver(x,y,u1(:,:,k)',u2(:,:,k)')
    title('Plots of the control which is actually the velocity');
    F(k) = getframe;
    writeVideo(vid,F(k));
end

close(vid);

% PDF

vid1 = VideoWriter('circle_pdf.avi');
vid1.FrameRate = 5;
open(vid1);

figure(2)
F1(N-1) = struct('cdata',[],'colormap',[]);

sol = fok_pl_adi2(u1,u2,rho,1);
for k = 1:N
    surf(x,y,sol(:,:,k)');
    view(0,90);
    colorbar
    title('PDF')
    F1(k) = getframe;
    writeVideo(vid1,F1(k));
end

close(vid1);