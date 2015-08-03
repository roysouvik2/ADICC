function val = V1(x,y,t1,t2)

% Formula for one trajectory.
val = ((x-t1(1))^2+(y-t1(2))^2);

% Formula for two trajectories.
%val = ((x-t1(1))^2+(y-t1(2))^2)*((x-t2(1))^2+(y-t2(2))^2);