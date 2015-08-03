function val = plot_traj()



for kt = 1:1
    [a,b, m, h, x, y, N, t, dt] = parameters(kt);
    for k = 1:N
        tr1 = traj1(t(k));
        %tr2 = traj2(t(k));
        figure(2)
        plot(tr1(1),tr1(2),'b');
        %pause(1)
        hold on;
        %plot(tr2(1),tr2(2),'b-');
        %pause(1)
        %hold on;
        axis([a b a b])
    end
end
        
    

