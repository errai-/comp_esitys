function [ accelerations ] = acceleration_total( locations,velocities,h)
% i) Finds a reasonable "important neighbors" listing
% ii) Calculates densities
% iii) Calculates finally accelerations, dumps all other temporary data

%% Initialization
m = 1;
particle_count = size(velocities,1);

neighbors = cell(particle_count,1);

densities = zeros(particle_count,1);
accelerations = zeros(size(velocities));
pressures = zeros(particle_count,1);

%% Naive neighbor listing, takes into account overlap
% TODO: intelligent neighbor listing

for i=1:particle_count
    neighbors{i} = i:particle_count;
end

%% Calculate densities

for i = 1:particle_count
    for j=neighbors{i},
        r = calc_distance_2D(locations(j,:),locations(i,:));
        tmpDens = m*cubic_spline_kernel(r,h);
        densities(i) = densities(i) + tmpDens;
        densities(j) = densities(j) + tmpDens;
    end
end
display(densities);

%% Calculate pressure
% using Eq.3.38 and 3.91 from the arxiv article, for gamma the value was
% taken from the book (usually used in HE simulations)
gamma = 3;
kappa = 1; % this entropic function is constant for isentropic flow (assumed)
for i = 1:particle_count
    internal_energy = kappa/(gamma - 1)*densities(i)^(gamma - 1);
    pressures(i) = (gamma - 1)*densities(i)*internal_energy;
end

%% Calculate accelerations:
% Loop over particles
for i=1:particle_count
    % Loop over neighbors not yet handled
    for j=neighbors{i},
        if (j == i), continue; end
        r = calc_distance_2D(locations(j,:),locations(i,:));
        direction = (-locations(j,:)+locations(i,:))/r;
        
        % calculate artificial viscosity (the uppercase pi_ij term)
        viscosity_condition = dot(velocities(i,:) - velocities(j,:),...
            locations(i,:) - locations(j,:)); % v_ij*r_ij
        if viscosity_condition < 0
            ave_density = (densities(i) + densities(j))/2;
            alpha = 1; % from the book
            beta = 10; % large value to prevent unphysical penetration in the exposion
            epsilon = 0.01; % from the arxiv article;
            mu = h*viscosity_condition/(calc_distance_2D(locations(i,:),...
                locations(j,:))^2 + epsilon*h^2);
            c_i = sqrt(gamma*kappa*densities(i)^(gamma - 1)); % sound speed dP/d(rho)
            c_j = sqrt(gamma*kappa*densities(j)^(gamma - 1));
            ave_c = (c_i + c_j)/2; % average sound speed
            viscosity = (-alpha*ave_c*mu + beta*mu^2)/ave_density; % standard SPH viscous term
        else
            viscosity = 0;
        end
        
        tmpAccel = m*direction*(pressures(i)/densities(i)^2+pressures(j)/densities(j)^2+...
        viscosity)*...
            cubic_spline_kernel_gradient(r,h);
        accelerations(i,:) = accelerations(i,:) - tmpAccel;
        accelerations(j,:) = accelerations(j,:) + tmpAccel;
    end
end

end

