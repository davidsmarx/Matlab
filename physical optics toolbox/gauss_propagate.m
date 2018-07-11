function [f, x, y] = gauss_propagate(z,wx,wy,lam,dx,dy,Nx,Ny)
% f = gauss_propagate(z,wx,wy,lam,dx,dy,Nx,Ny)
% inputs:
%    z = scalar distance from beam waist (zero phase)
%    wx, wy = gaussian beam radius (1/e^2 intensity) at beam waist in x and y
%    lam = wavelength in medium = vacuum wavelength / n
%    dx, dy = desired grid spacing for result field, same units as lam
%    Nx, Ny (optional) = # of grid points for output (default = 256 x 256)
% outputs:
%    f = matrix of complex field values
%    x, y = vectors of grid points for f
% method:
%    plane wave propagates the angular spectrum of a gaussian beam, then inverse
%    fft to get field

if nargin == 6,
   Nx = 256; Ny = 256;
elseif nargin ~= 8,
   error('usage: [f, x, y] = gauss_propagate(z,wx,wy,lam,dx,dy,Nx,Ny)');
end

du = lam/Nx/dx; dv = lam/Ny/dy;

u = [-Nx/2:Nx/2-1]'*du; v = [-Ny/2:Ny/2-1]'*dv;
[U, V] = meshgrid(u,v);
U2 = U.^2; clear U;
V2 = V.^2; clear V;

F = (-(pi/lam).^2)*(wx*wx*U2 + wy*wy*V2) + j*(2*pi.*z/lam)*sqrt(1-U2-V2);
clear U2 V2;
F = exp(F);

f = (pi*lam*lam*sqrt(wx*wy)) * fftshift(fft2(fftshift(F))) ./ (Nx*Ny);

x = [-Nx/2:Nx/2-1]'*dx;
y = [-Ny/2:Ny/2-1]'*dy;
return