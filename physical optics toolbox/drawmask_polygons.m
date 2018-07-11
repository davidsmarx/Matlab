function hh = drawmask_polygons(haxes,vertices,varargin)
% hh = drawmask_polygons(haxes,vertices,options)
%
% vertices is a 2 x n matrix. vertices(1,:) is x coordinate, vertices(2,:)
% is y coordinate.
%

options = varargin;

axes(haxes);

% lines, add first point to end to close polygon
xx = [vertices(1,:) vertices(1,1)];
yy = [vertices(2,:) vertices(2,1)];

% draw
hh(1) = line(xx,yy,options{:});
hh(2) = line(-xx,-yy,options{:});

