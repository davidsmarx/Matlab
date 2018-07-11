function hp = surfnogrid(x, y, z, w)
% hp = surfnogrid(x, y, z)
% hp = surfnogrid(x, y, z, w)
%
% x, y, z, w are vectors of the same length
% w determines the color (default w = z)
% 
% x, y do not need to be on a regular grid
% uses tri = delaunay(x,y) to create patches
%
% return value hp = array of patch handles

if ~exist('w','var') || isempty(w),
    w = z;
end

tri = delaunay(x, y);

hp = patch(x(tri)', y(tri)', z(tri)', w(tri)');
set(hp,'linestyle','none');
