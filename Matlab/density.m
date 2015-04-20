function [density] = density(particle_matrix,particle_index,h)
% Density at particle j. Based on Eq.3.39 of the arxiv paper

% particle matrix contains the neighbouring particle coordinates and has dimensions
% particle_count*dimensions where the dimension can be 2D. The matrix
% includes the particle j.

% particle_index is the index of the particle j in the matrix (place in the
% column)
% h is the smoothing length
% returns rho_j
% the particle mass is set to 1
m = 1;
density = 0;

particle_count = size(particle_matrix,1); 

for i = 1:particle_count % loop the particles
   r = calc_distance_2D(particle_matrix(particle_index,:),particle_matrix(i,:)); % the distance between the particles
   density = density + m*cubic_spline_kernel(r,h);
end


end