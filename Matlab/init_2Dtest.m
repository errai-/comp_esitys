clear all
N = 20;

dim = 2;
dim_sizes = [[-2 2];[-3 3]];
cell_size = 0.5;

t = 0;
m = 1;


[r, v, grid] = Initializing(N, dim, dim_sizes, cell_size);

r = random_box_2D(N, r, dim_sizes);

%eli meillä on alkusijainnit hiduille, ja kaikkien v = 0
%r että v ovat 2D matriiseja, x- ja y-komponenttijako on 2. indeksinä

%grid on cell, johon voi myöhemmin laittaa listoja hitu-indekseistä jotka
%on laatikoissa. Kirjotan myöhemmin funktiot 2D ja 3D
%jotka syö paikkakoordinaatit ja palauttaa, missä 
%laatikossa se sijainti on.
