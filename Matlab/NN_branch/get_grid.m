function [ grid_index, x_index, y_index ] = get_grid( location, dim_sizes, h )
%get_grid: returns the current location's grid index for 2D case as well as
%the "box indexes" in x and y directions

%h is grid size
%dim_sizes is matrix of corner points for the grid
%location is [x y]


x_len = (dim_sizes(1,2)-dim_sizes(1,1))/h;

x_index = ceil((dim_sizes(1,2) + location(1))/h);
y_index = ceil((dim_sizes(2,2) + location(2))/h);

%y_index = x_len*(ceil((dim_sizes(2,2) + location(2))/h)-1);

grid_index = x_index + x_len*(y_index-1);

end