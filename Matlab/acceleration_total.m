function [accelerations,hVals] = acceleration_total( locations,velocities,...
    hVals,m,kappa,gamma,rho_const)
% i) Finds a reasonable "important neighbors" listing
% ii) Calculates densities
% iii) Calculates finally accelerations, dumps all other temporary data

%% Initialization
particle_count = size(velocities,1);

% Pre-calculate the spline functions with neighbors
splines = cell(particle_count,1);
spline_gradients = cell(particle_count,1);
neighbors = cell(particle_count,1);

densities = zeros(particle_count,1);
accelerations = zeros(size(velocities));
pressure_per_rhos = zeros(particle_count,1);

%% Naive neighbor listing, takes into account overlap
% TODO: intelligent neighbor listing

for i=1:particle_count
    splineTmp = []; splineGradTmp = []; nbTmp = [];
    for j=i+1:particle_count
        rij = calc_distance_2D(locations(j,:),locations(i,:));
        hij = (hVals(i)+hVals(j))/2;
        if (rij < 2*hij),
            splineTmp = [splineTmp; cubic_spline_kernel(rij,hij)];
            splineGradTmp = [splineGradTmp; cubic_spline_kernel_gradient(rij,hij)*(locations(j,:)-locations(i,:))/rij];
            nbTmp = [nbTmp; j];
        end
    end
    splines{i} = splineTmp;
    spline_gradients{i} = splineGradTmp;
    neighbors{i} = nbTmp;
end

%% Calculate densities

for i = 1:particle_count
    densities(i) = densities(i) + m*cubic_spline_kernel(0,hVals(i)); 
    for j=1:size(neighbors{i},1),
        tmpDens = m*splines{i}(j);
        densities(i) = densities(i) + tmpDens;
        densities(neighbors{i}(j)) = densities(neighbors{i}(j)) + tmpDens;
    end
end

%% Calculate pressure
% using Eq.3.38 and 3.91 from the arxiv article, for gamma the value was
% taken from the book (usually used in HE simulations)
% The entropic function kappa is constant by assumption
for i = 1:particle_count
    pressure_per_rhos(i) = kappa*densities(i)^(gamma-1);
end

%% Calculate accelerations:
% Loop over particles
for i=1:particle_count
    % Loop over neighbors not yet handled
    for j=1:size(neighbors{i},1),
        
        % calculate artificial viscosity (the uppercase pi_ij term)
        viscosity_condition = dot(velocities(i,:) - velocities(j,:),...
            locations(i,:) - locations(j,:)); % v_ij*r_ij
        hij = (hVals(i)+hVals(j))/2;
        if viscosity_condition < 0
            ave_density = (densities(i) + densities(j))/2;
            alpha = 1; % from the book
            beta = 10; % large value to prevent unphysical penetration in the exposion
            epsilon = 0.01; % from the arxiv article;
            mu = hij*viscosity_condition/(calc_distance_2D(locations(i,:),...
                locations(j,:))^2 + epsilon*hij^2);
            c_i = sqrt(gamma*kappa*densities(i)^(gamma - 1)); % sound speed dP/d(rho)
            c_j = sqrt(gamma*kappa*densities(j)^(gamma - 1));
            ave_c = (c_i + c_j)/2; % average sound speed
            viscosity = (-alpha*ave_c*mu + beta*mu^2)/ave_density; % standard SPH viscous term
        else
            viscosity = 0;
        end
        
        tmpAccel = m*(pressure_per_rhos(i)/densities(i)+ ...
            pressure_per_rhos(j)/densities(j)+viscosity)* ...
            spline_gradients{i}(j,:);
        accelerations(i,:) = accelerations(i,:) + tmpAccel;
        accelerations(neighbors{i}(j),:) = accelerations(neighbors{i}(j),:) - tmpAccel;
    end

end
%% Update hVals

%deltas = zeros(particle_count,1);
%for i=1:particle_count
    %for j = 1:size(neighbors{i},1),
    %    delta = spline_gradients{i}(j,:)*(velocities(i,:)-velocities(j,:))';
    %    deltas(i) = deltas(i) + delta;
    %    deltas(j) = deltas(j) + delta;
    %end
    %hVals(i) = rho_const/sqrt(densities(i));%( 1 - deltas(i)*m/(2*densities(i)) );
%end
%hVals = (hVals+rho_const./sqrt(densities))/2;
display(mean(densities));

end
