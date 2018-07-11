function hh = drawmask_rects(haxes,hx,hy,xo,varargin)
% hh = drawmask_rects(haxes,hx,hy,xo,direction,options)
% hh = drawmask_rects(haxes,hx,hy,xo,rotangle,direction,options)
%
% direction = 'ns','ew','both'(default)

[rotangle, direction, options] = ParseInputs(varargin{:});

axes(haxes);

% vertices of rectangles
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

% draw
hh = [];
if strmatch(direction,{'ns','n2','both'}),
    hh = [hh line(yy,-xx,options{:})];
    hh = [hh line(-yy,xx,options{:})];
end
if strmatch(direction,{'ew','n1','both'}),
    hh = [hh line(xx,yy,options{:})];
    hh = [hh line(-xx,-yy,options{:})];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rotangle, direction, options] = ParseInputs(varargin)

if isempty(varargin),
    rotangle = 0;
    direction = 'both';
    options = {};
    return
end

switch class(varargin{1})
    case 'double',
        rotangle = varargin{1};
        if length(varargin) >= 2,
            direction = varargin{2};
            options = {varargin{3:end}};
        else,
            direction = 'both';
            options = {};
        end
    case 'char',
        rotangle = 0;
        direction = varargin{1};
        options = {varargin{2:end}};
end