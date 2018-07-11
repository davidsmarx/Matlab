function hh = drawmask_squareXrounded(haxes,fsmask,varargin)
% hh = drawmask_squareXrounded(haxes,fsmask,options)
%
% fsmask description:
%   fsmask(1) = length of each side for the ns square
%   fsmask(2) = radius from origin to center of squares for ew
%   fsmask(3) = length of each side for ew squares
%   fsmask(4) = (optional) radius of curvature for rounded corner of ew squares

axes(haxes);

% validate inputs
switch length(fsmask)
    case 3,
        le = fsmask(3);
        xc = fsmask(2)/sqrt(2);
        yc = fsmask(2)/sqrt(2);
        xx = [xc+le/2 xc+le/2 xc-le/2 xc-le/2 xc+le/2];
        yy = [yc+le/2 yc-le/2 yc-le/2 yc+le/2 yc+le/2];
        
    case 4,
        % define centers and points for 1st quadrant
        le = fsmask(3);
        rr = fsmask(4);
        xc = fsmask(2)/sqrt(2);
        yc = fsmask(2)/sqrt(2);
        xr = xc - le/2 + rr;
        yr = yc - le/2 + rr;
        t  = linspace(pi,3*pi/2,25);
        
        xx = [xr xc+le/2 xc+le/2 xc-le/2 xc-le/2 xr+rr*cos(t)];
        yy = [yc-le/2 yc-le/2 yc+le/2 yc+le/2 yr yr+rr*sin(t)];

    otherwise,
        error('field separator mask invalid parameters');
end

hh(1) = line( xx, yy, varargin{:});
hh(2) = line(-xx, yy, varargin{:});
hh(3) = line( xx,-yy, varargin{:});
hh(4) = line(-xx,-yy, varargin{:});
