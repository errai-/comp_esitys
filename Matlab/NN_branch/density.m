function [ densities ] = density( hVals,splines,neighbors, particle_count,m )
densities = zeros(particle_count,1);

for i = 1:particle_count
    densities(i) = densities(i) + m*cubic_spline_kernel(0,hVals(i)); 
    for j=1:size(neighbors{i},2),
        tmpDens = m*splines{i}(j);
        densities(i) = densities(i) + tmpDens;
        densities(neighbors{i}(j)) = densities(neighbors{i}(j)) + tmpDens;
    end
end

end

