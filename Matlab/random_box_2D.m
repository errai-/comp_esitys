function [ r ] = random_box_2D( N, r, dim_sizes )
%forms the starting positions for particles in 2D in simple
%manner

x_len = dim_sizes(1,2)-dim_sizes(1,1);
y_len = dim_sizes(2,2)-dim_sizes(2,1);


for i = 1:1:N
   r(i,1) = x_len*0.4 + x_len*0.2*rand() + dim_sizes(1,1);
   r(i,2) = y_len*0.4 + y_len*0.2*rand() + dim_sizes(2,1);
end
end


