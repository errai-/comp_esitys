function [neighbors,splines,spline_gradients] = ...
    neighbors_splinesNN(locations,hVals,particle_count, dim_sizes, grid, h_grid)

splines = cell(particle_count,1);
spline_gradients = cell(particle_count,1);
neighbors = cell(particle_count,1);

for i = 1:1:particle_count
    
    range = 2;
    %this is the search range, it *could* be changed for each particle
    %could be compared to the rij < 2*hij -condition
       
    grid_indexes = NNs(locations(i,:), dim_sizes, h_grid, range);
    idx = size(grid_indexes,1);
    grid_indexes = grid(grid_indexes);
    
    neighbors{i} = [grid_indexes{1:1:idx}];
    neighbors{i} = neighbors{i}(find(neighbors{i}>i));
    
end %endof for


for i = 1:1:particle_count
    distances = calc_distance_v2D(locations(neighbors{i},:),locations(i,:));
    h_vector = (hVals(i) + hVals(neighbors{i}))/2;
    splines{i} = zeros(1,size(neighbors{i},2));
    spline_gradients{i} = zeros(size(neighbors{i},2),2); %2 is hardcoded dimension!!!
    
    for j = 1:1:size(neighbors{i},2)
        splines{i}(j) = cubic_spline_kernel(distances(j),h_vector(j));
        
        %i
        %j
        %distances(j)
        %spline_gradients{i}(j)
        %cubic_spline_kernel_gradient(distances(j),h_vector(j))
        %locations(neighbors{i}(j))
        %cubic_spline_kernel_gradient(distances(j),h_vector(j))*(locations(neighbors{i}(j),:) - locations(i,:))/distances(j)
        
        spline_gradients{i}(j,:) = cubic_spline_kernel_gradient(distances(j),h_vector(j))*(locations(neighbors{i}(j),:) - locations(i,:))/distances(j);

    end %endof for
end %endof for


end

