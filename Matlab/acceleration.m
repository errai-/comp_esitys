function [ acceleration ] = acceleration( locations,velocities,densities,h )

particle_count = size(velocities,1);

for i=1:particle_count
    acceleration = [0,0];
end

end

