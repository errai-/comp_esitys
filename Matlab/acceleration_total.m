function accelerations = acceleration_total( locations,velocities,h)
% i) Finds a reasonable "important neighbors" listing
% ii) Calculates densities
% iii) Calculates finally accelerations, dumps all other temporary data

%% Initialization
m = 1;
particle_count = size(velocities,1);

rVals = cell(particle_count,1);
neighbors = cell(particle_count,1);

densities = zeros(particle_count,1);
accelerations = zeros(size(velocities));
pressures = zeros(particle_count,1);

%% Naive neighbor listing, takes into account overlap
% TODO: intelligent neighbor listing

for i=1:particle_count
    rTmp = []; nbTmp = [];
    for j=i+1:particle_count
        r = calc_distance_2D(locations(j,:),locations(i,:));
        if (r < 2*h),
            rTmp = [rTmp, r];
            nbTmp = [nbTmp, j];
        end
        rVals{i} = rTmp;
        neighbors{i} = nbTmp;
    end
end

%% Calculate densities

for i = 1:particle_count
    densities(i) = densities(i) + m*cubic_spline_kernel(0,h); 
    rTmp = rVals{i};
    rIdx = 0;
    for j=neighbors{i},
        rIdx = rIdx+1;
        tmpDens = m*cubic_spline_kernel(rTmp(rIdx),h);
        densities(i) = densities(i) + tmpDens;
        densities(j) = densities(j) + tmpDens;
    end
end

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
    rTmp = rVals{i};
    rIdx = 0;
    % Loop over neighbors not yet handled
    for j=neighbors{i},
        rIdx = rIdx + 1;
        direction = (-locations(j,:)+locations(i,:))/rTmp(rIdx);
        
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
        display(tmpAccel);
        display(pressures(i));
        display(densities(i));
        accelerations(i,:) = accelerations(i,:) - tmpAccel;
        accelerations(j,:) = accelerations(j,:) + tmpAccel;
    end
end


end
