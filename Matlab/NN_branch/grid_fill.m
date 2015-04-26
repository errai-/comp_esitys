function [ new_grid ] = grid_fill( locations, gsize, dim_sizes, h )
%GRID_UPDATE writes the grid listing anew completely

new_grid = cell(gsize,1);
N = size(locations,1);

for i = 1:1:N
   grid_idx = get_grid(locations(i,:),dim_sizes, h);
   new_grid{grid_idx} = [new_grid{grid_idx} i];      
end
   
    
    
end

