function A1 = angspec_propagate(A0,x,y,k,z)
% A1 = angspec_propagate(A0,x,y,k,z)

[Ny, Nx] = size(A0);

dx = mean(diff(x));
dy = mean(diff(y));

du = 1./dx./Nx;
dv = 1./dy./Ny;
u = [-Nx/2:Nx/2-1]'*du;
v = [-Ny/2:Ny/2-1]'*dv;
[U,V] = meshgrid(u,v);
ep = exp(j*sqrt(k.^2 - U.^2 - V.^2)*z);

At = ep.*fftshift(fft2(fftshift(A0)));

A1 = fftshift(ifft2(fftshift(At)));

return
