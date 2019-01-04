function com = calcCOM(im, varargin)
% com = calcCOM(im)
%
% simple function to calculate center of mass
% usuall input image is a bMask;
% default coordinates are 1:nc, 1:nr (i.e. 1-offset)

[nr, nc] = size(im);

x = CheckOption('x', (1:nc)', varargin{:});
y = CheckOption('y', (1:nr)', varargin{:});

[X, Y] = meshgrid(x,y);

com = sum([X(:).*im(:) Y(:).*im(:)])./sum(im(:));

