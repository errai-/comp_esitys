function [W] = cubic_spline_kernel_gradient(r,h)
% Cubic spline kerner derivative between two particles. Based on the arxiv paper pg.11

% r is abs(x_i - x_j) i.e. the distance between the particles
% h is the smoothing length
% returns grad(W_ij)
if r < 0, display('Misuse of kernel function'); end
x = r/h;

W = 1/(pi*h^3);

% This is also 3D
%W = 1/(pi*h^4);

%if x >= 0 && x <= 1
%    W = W*(9/4*x^2 - 3*x);
%end
%if x >= 1 && x <= 2
%    W = W*(-3/4)*(2 - x)^2;
%end

if x <= 1
    W = W*(3*x^2 - 4*x);
elseif x <= 2
    W = W*(-1)*(2 - x)^2;
else
    W = 0;
end

end