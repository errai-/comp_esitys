function [ accelerations ] = acceleration_total( locations,velocities,h )

accelerations = zeros(size(velocities));
particle_count = size(velocities,1);

densities = density_total(locations,h);

for i=1:particle_count
    accelerations(i,:) = acceleration( locations, velocities, densities, h );
end

end

