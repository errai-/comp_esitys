function [ locations, velocities ] = update_particles( locations, velocities, h, tStep )

% Calculate accelerations on the current step
accelerations = acceleration_total(locations,velocities,h);

% Half timestep values
temp_velocities = velocities+0.5*tStep*accelerations;
temp_locations = locations+0.5*tStep*velocities;
accelerations = acceleration_total(temp_locations,temp_velocities,h);

% Final step values
temp_velocities = velocities;
velocities = velocities+ tStep*accelerations;
locations = locations + 0.5*tStep*(temp_velocities+velocities);

end

