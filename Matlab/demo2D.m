clear all;
close all;

figure('Position', [1100, 100, 1049, 895]);

writerObj = VideoWriter('2D');
open(writerObj);

tInterval = 2;
m = 4;
kappa = 1; % A constant coefficient for pressures
gamma = 3; % Exponent in the pressure function
windowScale = 50;
N = 300;
len = 2;

%[locations,velocities] = colliding_walls(300,10);
locations = simple_random_2D( N, [len,len] );
%locations = simple_random_2D( N, [len/2,len/2] );
%locations = simple_random_2D( N, [len/4,4*len] );
velocities = zeros(size(locations));

% Approximate the initial density. This is used for initial h values
rho_init = N/( 1.2*len^2 );
hVals = ones( N, 1 )*(1/sqrt(4*pi))*sqrt( N*m / rho_init );

% Pre-calculate the spline functions with neighbors
[neighbors,splines,spline_gradients] = neighbors_splines(locations,hVals,N);
densities = density(hVals,splines,neighbors,N,m);
hVals = (1/sqrt(4*pi))*sqrt( N*m )./sqrt(densities);
hConst = (sqrt( densities )'*hVals)/N;

iter = 0;
printLen = 4;
tCurr = 0;
tStep = 0;
while tCurr < tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        display(tCurr);
        plot(locations(:,1),locations(:,2),'.');
        axis([-windowScale,windowScale,-windowScale,windowScale]);
        drawnow;
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    iter = iter+1;
    
    [locations,velocities,hVals,tStep] = update_particles(locations, ...
        velocities,hVals,m,kappa,gamma,hConst,tStep);
    tCurr = tCurr + tStep;
end

close(writerObj);
