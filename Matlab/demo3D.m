clear all;
close all;

figure('Position', [1100, 100, 1049, 895]);

writerObj = VideoWriter('3D');
open(writerObj);

tInterval = 5;
m = 1;
kappa = 1; % A constant coefficient for pressures
gamma = 5/3; % Exponent in the pressure function
windowScale = 50;
N = 300;
len = 1;

locations = simple_random_3D( N, [len,len,len] );
velocities = zeros(size(locations));

% Approximate the initial density. This is used for initial h values
rho_init = N/( 1.2*len^3 );
hVals = ones( N, 1 )*(1/sqrt(4*pi))*sqrt( N*m / rho_init );

% Pre-calculate the spline functions with neighbors
[neighbors,splines,spline_gradients] = neighbors_splines_3D(locations,hVals,N);
densities = density(hVals,splines,neighbors,N,m);
hVals = nthroot( (N*m*3/(32*pi))./densities, 3 );
hConst = (nthroot( densities, 3 )'*hVals)/N;

iter = 0;
printLen = 4;
tCurr = 0;
tStep = 0;
while tCurr < tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        
        display(tCurr);
        scatter3(locations(:,1),locations(:,2),locations(:,3),'.');
        axis([-windowScale,windowScale,-windowScale,windowScale,-windowScale,windowScale]);
        drawnow;
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    iter = iter+1;
    
    [locations,velocities,hVals,tStep] = update_particles_3D(locations, ...
        velocities,hVals,m,kappa,gamma,hConst,tStep);
    tCurr = tCurr + tStep;
end

close(writerObj);
