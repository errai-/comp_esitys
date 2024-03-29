function [W] = cubic_spline_kernel(r,h)
% Cubic spline kerner between two particles. Based on the arxiv paper pg.11

% r is abs(x_i - x_j) i.e. the distance between the particles
% h is the smoothing length
% returns the kernel W_ij
if r < 0, display('Misuse of kernel function'); end
x = r/h;

% NOTE: the original implementation was given in 3D. Saving it for
% a rainy day here in comments
%W = 1/(pi*h^3);
%if x >= 0 && x <= 1
%    W = W*(1 - 3/2*x^2 + 3/4*x^3);
%end
%if x >= 1 && x <= 2
%    W = W*1/4*(2 - x)^3;
%end

W = 1/(pi*h^2);

if  x <= 1
    W = W*(4/3 - 2*x^2 + x^3);
elseif x <= 2
    W = W*(1/3)*(2 - x)^3;
else    
    W = 0;
end

end

