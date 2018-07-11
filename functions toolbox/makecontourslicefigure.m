function h_fig = makecontourslicefigure(V,pixel_size,zz,sx,sy,sz)
% h_fig = makecontourslicefigure(V,pixel_size,zz,sx,sy,sz)
%
% V(:,:,z) is a grayscale image at focus = z
% pixel_size is the spacing between samples in the x-y plane, assumes dx=dy
% zz is a vector of positions in z (focus) where the images were recorded
% sx,sy,sz describe where to make the contourslices, see help contourslice

[ny,nx,nz] = size(V);
if length(zz) ~= nz, error('length(zz) ~= nz'); end

[X,Y,Z] = meshgrid([-(nx-1)/2:(nx-1)/2]*pixel_size,[-(ny-1)/2:(ny-1)/2]*pixel_size,zz);

h_fig = figure;
contourslice(X,Y,Z,V,sx,sy,sz)
grid,
set(gca,'view',[-37.5 30]), 
xlabel('X [\mum]'), ylabel('Y [\mum]'), zlabel('Depth [\mum]')
set(gca,'dataaspectratio',[1 1 1]),

set(gca,'xlim',[floor(X(1)) ceil(X(end))]), 
set(gca,'ylim',[floor(Y(1)) ceil(Y(end))])
