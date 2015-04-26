function [ locations, velocities, hVals, tStep, grid] = update_particles( locations, velocities, hVals, tStep,m,kappa,gamma,hConst, dim_sizes, grid, h_grid )

gsize = size(grid,1);
% Calculate accelerations on the current step
[accelerations,temp_hVals,tStep] = acceleration_total(locations, ...
    velocities,hVals,m,kappa,gamma,hConst, dim_sizes, grid, h_grid);

% Half timestep values
temp_velocities = velocities+0.5*tStep*accelerations;
temp_locations = locations+0.5*tStep*velocities;

temp_grid = grid_fill(temp_locations, gsize, dim_sizes, h_grid);


[accelerations,hVals] = acceleration_total(temp_locations, ...
    temp_velocities,temp_hVals,m,kappa,gamma,hConst, dim_sizes, temp_grid, h_grid);
%display(hVals);
% Final step values
temp_velocities = velocities;
velocities = velocities+ tStep*accelerations;
locations = locations + 0.5*tStep*(temp_velocities+velocities);

grid = grid_fill(locations, gsize, dim_sizes, h_grid);

end

