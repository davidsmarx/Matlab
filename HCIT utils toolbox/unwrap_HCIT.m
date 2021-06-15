function [phaunwrap, bMask] = unwrap_HCIT(pha, bMask, varargin)
% phaunwrap = unwrap_HCIT(pha, bMask)
%
% phase unwrapping specific to handle a CGI-like pupil with struts
% based on WFSC-CGI/Calibration/pr/util/unwrap.py
%

selem = CheckOption('selem', strel('disk',1), varargin{:}); % = [] to skip imerode(bMask)

if ~isempty(selem),
    bMask = imerode(bMask, selem);
end

% coordinate system with origin at center of mask
[com, x, y] = calcCOM(bMask);
x = x-com(1);
y = y-com(2);
[X, Y] = meshgrid(x,y); R = hypot(X, Y);
F = scatteredInterpolant(X(bMask), Y(bMask), pha(bMask), 'nearest', 'none');
rmax = max(R(bMask));
pha_interp = pha;
iinterp = R<rmax & ~bMask;
pha_interp(iinterp) = F(X(iinterp), Y(iinterp));
 
% call 2-d phase unwrap
phaunwrap = unwrap_phase(pha_interp);

% only bMask pixels are valid
phaunwrap(~bMask) = 0;



