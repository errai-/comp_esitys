function [distance] = calc_distance_2D(coordinates1,coordinates2)
% Calculates the distance between two particles in 2D case

distance_x = coordinates1(1) - coordinates2(1);
distance_y = coordinates1(2) - coordinates2(2);

distance = sqrt(distance_x^2 + distance_y^2);

end

