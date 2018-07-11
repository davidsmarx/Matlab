function [hh, xx, yy] = drawmask_fsmpoly(haxes,mask,bpoly,rotangle,direction,varargin)
% hh = drawmask_fsmpoly(haxes,mask,bpoly,rotangle,direction,options)
%
% mask(1) = distance to inside line
% mask(2) = dy
%
% bpoly = 2 x n, bpoly(1,:) = x-coordinates, bpoly(2,:) = y
% bpoly is the corner cube minimum aperture
%
% direction = 'ns','ew','both'(default)
%
% hh = array of handles to line objects
% xx, yy = vertices of one of the rectangles

axes(haxes);

% vertices of rectangles
XB = max(sqrt(sum(bpoly.^2)));
hx = XB - mask(1);
xo = (XB + mask(1))/2;
hy = mask(2);
x1 = xo - hx/2;
x2 = xo + hx/2;
y1 = -hy/2;
y2 =  hy/2;

% lines
xx = [x1 x2 x2 x1 x1];
yy = [y1 y1 y2 y2 y1];

% rotate
R = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
atmp = R * [xx; yy];
xx = atmp(1,:);
yy = atmp(2,:);

% draw rectangles
hh = [];
if strmatch(direction,{'ns','n2','both'}),
    hh = [hh line(yy,-xx,varargin{:})];
    hh = [hh line(-yy,xx,varargin{:})];
end
if strmatch(direction,{'ew','n1','both'}),
    hh = [hh line(xx,yy,varargin{:})];
    hh = [hh line(-xx,-yy,varargin{:})];
end

% draw corner cube aperture poly
bpolyc = [bpoly bpoly(:,1)]; % close polygon
hh = line(bpolyc(1,:),bpolyc(2,:),varargin{:});
hh = line(-bpolyc(1,:),-bpolyc(2,:),varargin{:});