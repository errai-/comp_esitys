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

%% Calculate accelerations:
% Loop over particles
for i=1:particle_count
    % Loop over neighbors not yet handled
    for j=neighbors{i},
        if (j == i), continue; end
        r = calc_distance_2D(locations(j,:),locations(i,:));
        direction = (-locations(j,:)+locations(i,:))/r;
        tmpAccel = m*direction*(1/densities(i)+1/densities(j))*...
            cubic_spline_kernel_gradient(r,h);
        accelerations(i,:) = accelerations(i,:) - tmpAccel;
        accelerations(j,:) = accelerations(j,:) + tmpAccel;
    end
end

end

