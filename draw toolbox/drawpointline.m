function [h, p1, polygon] = drawpointline(haxes,p0,dx,dy,varargin)
% h = drawpointline(haxes,p0,dx,dy,varargin)
%
% draw lines beginning at p0 = [x y]
% dx = list of x-lengths
% dy = list of y-lengths
%
% h = handle to line
% p1 = [x y] of last point
% polygon = [x(:) y(:)] list of vertices that make the polygon

axes(haxes);

xx = cumsum([p0(1); dx(:)]);
yy = cumsum([p0(2); dy(:)]);

h = line(xx,yy,varargin{:});

p1 = [xx(end) yy(end)];
polygon = [xx(:) yy(:)];