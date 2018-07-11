function h = drawrect(haxes,lenx,leny,varargin)
% h = drawrect(haxes,lenx,leny,x0,y0,options)
%
% haxes = handle of axes where to draw the circle
% lenx = length of rectangle in x-direction
% leny = length of rectangle in y-direction
% x0, y0 (optional) = center of circle
% options are any valid option used in line()
%
% output h = line handle

% validate optional input arguments and set defaults
[haxes,lenx,leny,x0,y0,optionlist] = ValidateInputs(haxes,lenx,leny,varargin{:});

axes(haxes);

x1 = x0-lenx/2; x2 = x0+lenx/2;
y1 = y0-leny/2; y2 = y0+leny/2;
x = [x1 x2 x2 x1 x1];
y = [y1 y1 y2 y2 y1];

h = line(x,y,optionlist{:});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [haxes,lenx,leny,x0, y0, optionlist] = ValidateInputs(ha,lx,ly,varargin)

% defaults
if ~isempty(ha), haxes = ha; else, haxes = gca; end
if lx < 0, error('drawrect: length < 0'); else, lenx = lx; end
if ly < 0, error('drawrect: length < 0'); else, leny = ly; end

x0 = 0; y0 = 0; optionlist = {};

narg = length(varargin);
if narg >= 1 & ~isempty(varargin{1}), x0 = varargin{1}; end
if narg >= 2 & ~isempty(varargin{2}), y0 = varargin{2}; end
if narg >= 3 & ~isempty(varargin{3}), optionlist = {varargin{3:end}}; end

return