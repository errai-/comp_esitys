function [locations,velocities] = colliding_walls(N,len)
    % collision of two asymmetric walls

    locations = zeros(N,2);
    locations_y = linspace(-len,len,N/2)';
    locations(1:N/2,2) = locations_y;
    locations(N/2+1:end,2) = locations_y + 5;
    locations(1:N/2,1) = -1;
    locations(N/2+1:end,1) = 1;
    
    velocities = zeros(size(locations));
    velocities(1:N/2,1) = 10;
    velocities(N/2+1:end,1) = -10;

end

