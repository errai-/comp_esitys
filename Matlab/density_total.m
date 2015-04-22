function [densities] = density_total(locations,h)

particle_count = size(locations,1); 
densities = zeros(particle_count,1);

% TODO: calculate nearest neighbors
for i = 1:particle_count
    densities(i) = density(locations,i,h);
end

end