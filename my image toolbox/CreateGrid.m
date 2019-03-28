function [x, y, X, Y, R, T] = CreateGrid(Ima, dx, dy, varargin)
% [x, y, X, Y, R, T] = CreateGrid(Ima, dx, dy)
% [x, y, X, Y, R, T] = CreateGrid(N, dx, dy)
% [x, y, X, Y, R, T] = CreateGrid([Ny, Nx], dx, dy)
% [x, y, X, Y, R, T] = CreateGrid([Ny, Nx], [dx, dy])
% dx (optional)
% dy (optional)
% [...] = CreateGrid(N, dx, dy, 'origin','lowerleft')

OriginLoc = CheckOption('origin','center',varargin{:});

switch nargin,
    case 1
        dx = 1;
        dy = 1;
    case 2,
        if length(dx) == 1,
            dy = dx;
        else
            dy = dx(2);
            dx = dx(1);
        end
end
        
if isscalar(Ima),
    ny = Ima; nx = Ima;
elseif isequal(size(Ima),[1 2]),
    ny = Ima(1); nx = Ima(2);
else
    [ny, nx] = size(Ima);
end

switch OriginLoc
    case 'center',
        x = dx * (ceil(-nx/2):ceil(nx/2-1))';
        y = dy * (ceil(-ny/2):ceil(ny/2-1))';
    case 'lowerleft',
        x = (0:nx-1)'*dx;
        y = (0:ny-1)'*dy;
    otherwise
        error(['unknown origin option: ' OriginLoc]);
end

[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
T = atan2(Y,X);