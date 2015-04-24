function [accelerations,hVals] = acceleration_total( locations,velocities,...
    hVals,m,kappa,gamma,tStep)
% i) Finds a reasonable "important neighbors" listing
% ii) Calculates densities
% iii) Calculates finally accelerations, dumps all other temporary data

%% Initialization
particle_count = size(velocities,1);

accelerations = zeros(size(velocities));
deltas = zeros(particle_count,1);
pressure_per_rhos = zeros(particle_count,1);

%% Naive neighbor listing, takes into account overlap
% TODO: intelligent neighbor listing
[neighbors,splines,spline_gradients] = neighbors_splines(locations,hVals,particle_count);

%% Calculate densities
densities = density(hVals,splines,neighbors,particle_count,m);

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
         % For hVal calculation
        delta = spline_gradients{i}(j,:)*(velocities(i,:)-velocities(j,:))';
        deltas(i) = deltas(i) - delta;
        deltas(j) = deltas(j) + delta;
        
        % calculate artificial viscosity (the uppercase pi_ij term)
        viscosity_condition = (velocities(i,:)-velocities(j,:))* ...
            (locations(i,:) - locations(j,:))'; % v_ij*r_ij
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
hVals = hVals.*( ones(particle_count,1) - tStep*deltas*m./(2*densities) );

end
