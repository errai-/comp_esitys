clear all;
figure('Position', [100, 100, 1049, 895]);

tStep = 0.01;
tInterval = 5;
m = 1;

% TODO: generate particles smoothly on some surface area

% An initial grid of locations and velocities
[X,Y] = meshgrid(0:0.1:1,0:0.1:1);
A = 1; % Area of the particle mesh

%locations = zeros(11*11,2);
locations = rand(11*11,2);
%locations(:,1) = reshape(X,[size(locations,1),1]);
%locations(:,2) = reshape(Y,[size(locations,1),1]);

velocities = zeros(size(locations));

% Approximates 10 particles in the neighborhood of a particle
particle_count = size(locations,1);
hInit = nthroot(10*A/(pi*particle_count),3)/2;

iter = 0;
printLen = round(0.1/tStep);

for tCurr = 0:tStep:tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        display(tCurr);
        plot(locations(:,1),locations(:,2),'*');
        axis([-2,2,-2,2]);
        pause;
    end
    iter = iter+1;
    
    [locations,velocities] = update_particles(locations,velocities,hInit,tStep);
end