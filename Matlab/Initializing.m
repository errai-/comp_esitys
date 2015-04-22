function [r, v, grid] = Initializing( N, dim, dim_sizes, cell_size )
%Alustaa hitujen sijainnit, nopeudet ja luo cell:n johon indeksoida missä
%laatikossa hidut ovat

%INPUT
%dim = montako dimensiota, esim 1, 2, 3
%dim_sizes = kokoa (dim, 2) esim [-1 1] [-2 2] -tyylinen laatikon kokolista
%cell_size = hitujen laskentacellien koko
%N = montako hitua

%OUTPUT
r = zeros(N, dim);
v = zeros(N, dim);

%grid on yksi pitkä cell-vektori, jonka pituus on laatikoiden määrä koko
%N-ulotteisessa avaruudessa

grid_size = 1;
for i = 1:1:dim
    grid_size = grid_size*(dim_sizes(i,2)-dim_sizes(i,1))/cell_size;
end
grid = cell(grid_size,1);

end


