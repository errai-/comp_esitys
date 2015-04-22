clear all;
figure(1);

tStep = 0.01;
tInterval = 5;

positions = [0,0;0,1;0,0.5;-0.5,0];
velocities = [0,0;1,0;0,0;0,0];

iter = 0;
printLen = round(0.01/tStep);

for tCurr = 0:tStep:tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        display(tCurr);
        plot(positions(:,1),positions(:,2),'*');
        pause;
    end
    iter = iter+1;
    
    [positions,velocities] = update_particles(positions,velocities,tStep);
end