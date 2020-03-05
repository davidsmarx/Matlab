function [com, x, y] = calcCOM(x, y, im, varargin)
% [com, x, y] = calcCOM(x, y, im)
% [com, x, y] = calcCOM([], [], im) (default (x, y) 1:nc, 1:nr
% [com, x, y] = calcCOM(im)
%
% simple function to calculate center of mass
% usual input image is a bMask;
% default coordinates are 1:nc, 1:nr (i.e. 1-offset)

if nargin == 1,
    im = x;
    x = [];
end
   
[nr, nc] = size(im);

if ~exist('x','var') || isempty(x),
    x = 1:nc;
end

if ~exist('y','var') || isempty(y),
    y = 1:nr;
end

if isscalar(x),
    x = x + (0:nc-1);
end
if isscalar(y),
    y = y + (0:nr-1);
end

% x = CheckOption('x', (1:nc)', varargin{:});
% y = CheckOption('y', (1:nr)', varargin{:});

[X, Y] = meshgrid(x,y);

com = sum([X(:).*im(:) Y(:).*im(:)])./sum(im(:));

