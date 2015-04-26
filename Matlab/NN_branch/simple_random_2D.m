function [ r ] = simple_random_2D( N, dim_sizes )
%forms the starting positions for particles in 2D in simple
%manner

r(:,1) = dim_sizes(1)*(rand(N,1)-0.5*ones(N,1));
r(:,2) = dim_sizes(2)*(rand(N,1)-0.5*ones(N,1));
