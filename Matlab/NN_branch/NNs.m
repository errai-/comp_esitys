function [ grid_indexes ] = NNs( location, dim_sizes, h, range )
%Nearest Neightbours at location, in grid, in 2D, including the particle
%itself

%range is in boxes

x_len = (dim_sizes(1,2)-dim_sizes(1,1))/h;
y_len = (dim_sizes(2,2)-dim_sizes(2,1))/h;

[~, x_c, y_c] = get_grid(location, dim_sizes, h);


grid_indexes = zeros((2*range+1)*(2*range+1),1);

idx = 0;
x_min = max([x_c-range,1]);
x_max = min([x_len, x_c+range]);

y_min = max([y_c-range,1]);
y_max = min([y_len, y_c+range]);

for r = x_min:1:x_max
        x_index = r;
        
        for s = y_min:1:y_max
                idx = idx +1;
                y_index = s;
                grid_indexes(idx) = x_index + x_len*(y_index-1);
            
        end %endof for
    
end %endof for

grid_indexes = grid_indexes(1:idx);


end







