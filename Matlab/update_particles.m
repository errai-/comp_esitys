function [ locations, velocities, hVals ] = update_particles( locations, velocities, hVals, tStep,m,kappa,gamma )

% Calculate accelerations on the current step
[accelerations,temp_hVals] = acceleration_total(locations, ...
    velocities,hVals,m,kappa,gamma);
plot(locations(:,1),locations(:,2),'.');
axis([-20,20,-20,20]);
pause;

% Half timestep values
temp_velocities = velocities+0.5*tStep*accelerations;
temp_locations = locations+0.5*tStep*velocities;
[accelerations,hVals] = acceleration_total(temp_locations, ...
    temp_velocities,temp_hVals,m,kappa,gamma);

% Final step values
temp_velocities = velocities;
velocities = velocities+ tStep*accelerations;
locations = locations + 0.5*tStep*(temp_velocities+velocities);

end

