function [x, y, X, Y, R, T] = CreateGrid(Ima, dx, dy, varargin)
% [x, y, X, Y, R, T] = CreateGrid(Ima, dx, dy)
% [x, y, X, Y, R, T] = CreateGrid(N, dx, dy)
% [x, y, X, Y, R, T] = CreateGrid([Ny, Nx], dx, dy)
% [x, y, X, Y, R, T] = CreateGrid([Ny, Nx], [dx, dy])
% dx (optional)
% dy (optional)
% options:
%    'origin', 'center' (default), '0-offset', '1-offset', 'center-halfpixel'
%    'xy_offset', [xo yo], after origin is set, x => x - xo, y=> y - yo

OriginLoc = CheckOption('origin','center',varargin{:});
xyo = CheckOption('xy_offset', [0 0], varargin{:});

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

if isempty(dx), dx = 1; end
if isempty(dy), dy = 1; end
        
if isscalar(Ima),
    ny = Ima; nx = Ima;
elseif isequal(size(Ima),[1 2]),
    ny = Ima(1); nx = Ima(2);
else
    sizeIma = size(Ima); % could be [nr nc] or [nr nc 3]
    ny = sizeIma(1);
    nx = sizeIma(2);
end

switch OriginLoc
    case 'center',
        x = dx * (ceil(-nx/2):ceil(nx/2-1))';
        y = dy * (ceil(-ny/2):ceil(ny/2-1))';
    case {'lowerleft', '0-offset'}
        x = (0:nx-1)'*dx;
        y = (0:ny-1)'*dy;
    case '1-offset'
        x = (1:nx)'*dx;
        y = (1:ny)'*dy;
    case 'center-halfpixel',
        % grid is half-pixel offset, so that origin is boundary between
        % pixels, if even # pixels
        % e.g. N=4, x = [-1.5 -0.5 0.5 1.5]
        x = dx * (0.5 + (ceil(-nx/2):ceil(nx/2-1))');
        y = dy * (0.5 + (ceil(-ny/2):ceil(ny/2-1))');
    otherwise
        error(['unknown origin option: ' OriginLoc]);
end

x = x - xyo(1);
y = y - xyo(2);

[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
T = atan2(Y,X);