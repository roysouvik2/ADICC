function val = plot_traj()

[a,b, m, h, x, y, N, t, dt] = parameters(1);


% rho = zeros(m+1,m+1);
% % Assigning initial pdf.
% for i = 1:m+1
%     for j = 1:m+1
%         %if(i==(m+1)/2 && j==(m+1)/2) 
%         %if((x(i)-1)^2 +(y(j)-1)^2 <1^2 )
%          rho(i,j) = exp((-(x(i))^2-(y(j))^2)/0.5);
%             %rho(i,j) = 1.0;
%         %else
% %             rho(i,j) = 0;
% %         end
%      end
% end
% rho = rho/(sum(sum(rho))*h^2);
% sol = fok_pl_adi2(u1,u2,rho,1);
% for k = 1:N
%     figure(1)
%     %mesh(sol(:,:,k)');
%     surf(x,y,sol(:,:,k)');
%     colorbar
%     view(0,90);
%     pause(1);
% end

% Plotting the gradient
% 
for k = 1:N-1
    figure(2)
    surf(x,y,u1(:,:,k)')
    colorbar;
    view(0,90)
    pause(1)
    figure(3)
    surf(x,y,u2(:,:,k)')
    colorbar;
    view(0,90)
    pause(1)
end