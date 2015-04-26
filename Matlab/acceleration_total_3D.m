function [accelerations,tStep] = acceleration_total_3D( locations,velocities,...
    hVals,m,kappa,gamma,tStep)
% i) Finds a reasonable "important neighbors" listing
% ii) Calculates densities
% iii) Calculates finally accelerations, dumps all other temporary data

%% Initialization
particle_count = size(velocities,1);

accelerations = zeros(size(velocities));
pressure_per_rhos = zeros(particle_count,1);
tStep_candidates = zeros(particle_count,1);

% viscosity constants
alpha = 1; % from the book
beta = 10; % large value to prevent unphysical penetration in the exposion
epsilon = 0.01; % from the arxiv article;

%% Naive neighbor listing, takes into account overlap
% TODO: intelligent neighbor listing
[neighbors,splines,spline_gradients] = neighbors_splines_3D(locations,hVals,particle_count);

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
        % calculate artificial viscosity (the uppercase pi_ij term)
        viscosity_condition = (velocities(i,:)-velocities(j,:))* ...
            (locations(i,:) - locations(j,:))'; % v_ij*r_ij
        hij = (hVals(i)+hVals(j))/2;
        if viscosity_condition < 0
            ave_density = (densities(i) + densities(j))/2;
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

%% Optimizing the time step with CFL criterion and force condition
% Using a global time step in order to avoid information lag between the
% particles. Eq.3.163 and 3.164 from arxiv paper (the generalized criteria
% is not valid when starting with zero velocities)

if (tStep == 0),
    for i = 1:particle_count
        % CFL criterion
        norm_grad_V = norm(accelerations(i,:));
        c_i = sqrt(gamma*kappa*densities(i)^(gamma - 1));
        if norm_grad_V < 0
            tStep_CFL = 0.3*hVals(i)/(c_i + hVals(i)*norm_grad_V + ...
                1.2*(alpha*c_i + beta*hVals(i)*norm_grad_V));
        else
            tStep_CFL = 0.3*hVals(i)/(c_i + hVals(i)*norm_grad_V);
        end
        % Force condition
        tStep_force = 0.3*sqrt(hVals(i)/norm(accelerations(i,:)));

        if 10*tStep_CFL < 0.1 % cheating (in the beginning the CFL time step is very small)
            tStep_CFL = 100*tStep_CFL; % change to 10 for precision
        end
        tStep_candidates(i) = min(tStep_CFL,tStep_force);
    end
    
    tStep = min(tStep_candidates); % the global time step

    % mitigating unwanted tStep behaviour close to the break down point
    if isfinite(tStep) == 0 % capture inf and NaN
       tStep = 0.01;
    elseif tStep > 0.1
       tStep = 0.001; 
    elseif tStep <= 0
       tStep = 0.001; 
    end
end
    
end
