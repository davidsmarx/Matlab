function h = drawcirc(haxes,r,varargin)
% h = drawcirc(haxes,r,x0,y0,options)
%
% haxes = handle of axes where to draw the circle, use [], or gca for
%     current axes
% r = radius of circle
% x0, y0 (optional) = center of circle
% options are any valid option used in line()
%
% output h = line handle

[haxes, r, x0, y0, optionlist] = ValidateInputs(haxes,r,varargin{:});

t = linspace(0,2*pi);

x = r.*cos(t)+x0;
y = r.*sin(t)+y0;

axes(haxes);

h = line(x,y,optionlist{:});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [haxes, r, x0, y0, optionlist] = ValidateInputs(htmp,rtmp,varargin)

% defaults
haxes = gca; r = 1; x0 = 0; y0 = 0; optionlist = {};

if ~isempty(htmp), haxes = htmp; end
if rtmp > 0, r = rtmp; end

if length(varargin) >= 1 & ~isempty(varargin{1}), x0 = varargin{1}; end
if length(varargin) >= 2 & ~isempty(varargin{2}), y0 = varargin{2}; end
if length(varargin) >= 3, optionlist = {varargin{3:end}}; end

return