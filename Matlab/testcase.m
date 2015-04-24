clear all;
close all;

figure('Position', [1100, 100, 1049, 895]);

tStep = 0.01;
tInterval = 100;
m = 1;
kappa = 1; % A constant coefficient for pressures
gamma = 3; % Exponent in the pressure function
windowScale = 30;
N = 300;
len = 3;

locations = simple_random_2D( N, [len,len] );
velocities = zeros(size(locations));

% Approximate the initial density. This is used for initial h values
rho_init = N/( 1.2*len^2 );
hVals = ones( N, 1 )*(1/sqrt(4*pi))*sqrt( N*m / rho_init );

% Pre-calculate the spline functions with neighbors
splines = cell(N,1);
neighbors = cell(N,1);
densities = zeros(N,1);

for i=1:N
    splineTmp = []; splineGradTmp = []; nbTmp = [];
    for j=i+1:N
        rij = calc_distance_2D(locations(j,:),locations(i,:));
        hij = (hVals(i)+hVals(j))/2;
        if (rij < 2*hij),
            splineTmp = [splineTmp; cubic_spline_kernel(rij,hij)];
            nbTmp = [nbTmp; j];
        end
    end
    splines{i} = splineTmp;
    spline_gradients{i} = splineGradTmp;
    neighbors{i} = nbTmp;
end
for i = 1:N
    densities(i) = densities(i) + m*cubic_spline_kernel(0,hVals(i)); 
    for j=1:size(neighbors{i},1),
        tmpDens = m*splines{i}(j);
        densities(i) = densities(i) + tmpDens;
        densities(neighbors{i}(j)) = densities(neighbors{i}(j)) + tmpDens;
    end
end
hVals = (1/sqrt(4*pi))*sqrt( N*m )./sqrt(densities);

iter = 0;
printLen = 1;
for tCurr = 0:tStep:tInterval
    % Loop logistics
    if ( mod(iter,printLen) == 0 )
        display(tCurr);
        plot(locations(:,1),locations(:,2),'.');
        axis([-windowScale,windowScale,-windowScale,windowScale]);
        drawnow;
    end
    iter = iter+1;
    
    [locations,velocities,hVals] = update_particles(locations, ...
        velocities,hVals,tStep,m,kappa,gamma);
end
