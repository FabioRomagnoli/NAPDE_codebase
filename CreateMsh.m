function [msh] = CreateMsh(x_int, x0, coordinates)
%CREATEMSH Summary of this function goes here
%   Detailed explanation goes here
dx_int=diff(x_int);
np = numel(dx_int); % number of points= number of cells
x=(x_int(2:end)+x_int(1:end-1))/2; %center of the cell

dx=diff(x);% distance among the cell centers
dx=[dx_int(1)/2;dx;dx_int(end)/2];% add distance with the ghost node

x = x + x0;
x_int = x_int + x0;

if strcmp (coordinates, "cartesian")
    msh.Vol = dx_int;
    msh.Area = ones(np+1, 1);
elseif strcmp (coordinates, "cylindrical")
    msh.Vol = dx_int .* (x_int(1:end-1) + x_int(2:end)) / 2;
    msh.Area = x_int;
end
msh.inv_Vol = 1 ./ msh.Vol;

msh.w = dx(1) / dx(2);
msh.e = dx(end) / dx(end-1);
msh.dx = dx;
msh.dx_int = dx_int;
msh.x = x;
msh.x_int = x_int;
msh.np = np;
msh.Area_over_dx = msh.Area ./ msh.dx;

end
