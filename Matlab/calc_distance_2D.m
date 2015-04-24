function [distance] = calc_distance_2D(coordinates1,coordinates2)
% Calculates the distance between two particles in 2D case

distance = sqrt( sum( (coordinates1-coordinates2).^2 ) );

end

