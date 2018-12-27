function [seed_coord]= voxel_coord_matrix_function(vol_dim, coord_spacing)
% Christopher Zahasky
% 12/20/2017

x = [1:coord_spacing:vol_dim];
y = [1:3:vol_dim];
z = [1:3:vol_dim];
[X, Y, Z] = meshgrid(x,y,z);
seed_coord = [X(:), Y(:), Z(:)];
