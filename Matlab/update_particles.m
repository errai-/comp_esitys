function [ locations, velocities, hVals, tStep] = update_particles( locations, velocities, hVals,m,kappa,gamma,hConst,tStep )

particle_count = size(velocities,1);

% Calculate accelerations on the current step
[accelerations,tStep] = acceleration_total(locations, ...
    velocities,hVals,m,kappa,gamma,tStep);

% Half timestep values
temp_velocities = velocities+0.5*tStep*accelerations;
temp_locations = locations+0.5*tStep*velocities;
accelerations = acceleration_total(temp_locations, ...
    temp_velocities,hVals,m,kappa,gamma,tStep);

% Final step values
temp_velocities = velocities;
velocities = velocities+ tStep*accelerations;
locations = locations + 0.5*tStep*(temp_velocities+velocities);

% Calculate hVals for the final timestep - Implicit euler or hard reset
[neighbors,splines] = neighbors_splines(locations,hVals,particle_count);
densities = density(hVals,splines,neighbors,particle_count,m);
hVals = hConst./sqrt( densities );

end

