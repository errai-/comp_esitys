function [distances] = calc_distance_2D(vector,coordinates)
% Calculates the distance between two particles in 2D case


distances_x = vector(:,1) - coordinates(1);
distances_y = vector(:,2) - coordinates(2);

distances = sqrt(distances_x.^2+distances_y.^2);

end
