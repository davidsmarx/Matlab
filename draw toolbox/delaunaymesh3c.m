function hout = delaunaymesh3c(x, y, z, c, varargin)
% h = delaunaymesh3c(x, y, z, c, optional properties)
% 
% c determines color of each path:
%   if c is same size as x and y, then each node is colored according to
%   the colormap with index c
%   if c is a string, e.g. 'r','g',..., then the whole plot is that color
%   if c is [r g b], then the whole plot is that color
% if no optional properties, then will use 'facecolor','interp'
%
% return value h = handle to the patch object

% default c = z
if nargin == 3,
    c = z;
end

% interpret the color argument
if isstr(c),
end

% set options according to varargin
if isempty(varargin)
    optionslist = {'facecolor','interp'};
else
    optionslist = varargin;
end

tri = delaunay(x, y);

% set patch colors based on c input
if length(c) == length(x),
    cpat = c(tri)';
else
    cpat = c;
end

h = patch(x(tri)',y(tri)',z(tri)',cpat,optionslist{:});

if nargout > 0,
    hout = h;
end

grid on
colorbar
view(3)