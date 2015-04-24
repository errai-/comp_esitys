function [neighbors,splines,spline_gradients] = ...
    neighbors_splines(locations,hVals,particle_count)
splines = cell(particle_count,1);
spline_gradients = cell(particle_count,1);
neighbors = cell(particle_count,1);

for i=1:particle_count
    splineTmp = []; splineGradTmp = []; nbTmp = [];
    for j=i+1:particle_count
        rij = calc_distance_2D(locations(j,:),locations(i,:));
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