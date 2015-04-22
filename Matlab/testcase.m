clear all;
close all;

figure('Position', [100, 100, 1049, 895]);

tStep = 0.01;
tInterval = 100;
m = 1;
kappa = 1; % A constant coefficient for pressures
gamma = 3; % Exponent in the pressure function
windowScale = 50;
N = 16;

hInit = 1;

locations = simple_random_2D( N, [2,2] );
velocities = zeros(size(locations));

iter = 0;
printLen = 100;
for tCurr = 0:tStep:tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        display(tCurr);
        plot(locations(:,1),locations(:,2),'.');
        axis([-windowScale,windowScale,-windowScale,windowScale]);
        pause;
    end
    iter = iter+1;
    
    [locations,velocities] = update_particles(locations,velocities,hInit,tStep,m,kappa,gamma);
end