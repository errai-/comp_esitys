clear all;
close all;

figure('Position', [1100, 100, 1049, 895]);

tStep = 0.005;
tInterval = 100;
m = 1;
kappa = 1; % A constant coefficient for pressures
gamma = 3; % Exponent in the pressure function
windowScale = 30;
N = 300;
len = 2;

%[locations,velocities] = colliding_walls(300,10);
locations = simple_random_2D( N, [len,len] );
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
printLen = 1;
tCurr = 0;
while tCurr < tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        %display(tCurr);
        plot(locations(:,1),locations(:,2),'.');
        axis([-windowScale,windowScale,-windowScale,windowScale]);
        drawnow;
    end
    iter = iter+1;
    
    % Hard reset for h values every 20 steps
    if (mod(iter,20)==0 && iter~=0),
        hTmp = hConst;
    else
        hTmp = 0;
    end
    [locations,velocities,hVals,tStep] = update_particles(locations, ...
        velocities,hVals,tStep,m,kappa,gamma,hTmp);
    tCurr = tCurr + tStep;
end
