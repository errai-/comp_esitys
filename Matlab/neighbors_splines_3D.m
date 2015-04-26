function [neighbors,splines,spline_gradients] = ...
    neighbors_splines_3D(locations,hVals,particle_count)
splines = cell(particle_count,1);
spline_gradients = cell(particle_count,1);
neighbors = cell(particle_count,1);

for i=1:particle_count
    splineTmp = []; splineGradTmp = []; nbTmp = [];
    distances = calc_distance_v3D(locations(i+1:end,:),locations(i,:));
    for j=i+1:particle_count
        rij = distances(j-i,:);
        hij = (hVals(i)+hVals(j))/2;
        if (rij < 2*hij),
            splineTmp = [splineTmp; cubic_spline_kernel(rij,hij)];
            splineGradTmp = [splineGradTmp; cubic_spline_kernel_gradient(rij,hij)*(locations(j,:)-locations(i,:))/rij];
            nbTmp = [nbTmp; j];
        end
    end
    splines{i} = splineTmp;
    spline_gradients{i} = splineGradTmp;
    neighbors{i} = nbTmp;
end

end