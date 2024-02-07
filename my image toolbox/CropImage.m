function [Imcrop, bMaskcrop] = CropImage(Im, bMask, xycent, xwid, ywid, varargin)
% [Imcrop, bMaskcrop] = CropImage(Im, bMask, xycent, xwid, ywid)
% simple square cropping of the image and its bMask
%
% xycent = [xc yc] in pixels, with [0, 0] at the center of the image array
% xwid, ywid are in pixels
%
% crop to array such that:
%   ixuse = x >= -xwid/2 & x < xwid/2;
%   iyuse = y >= -ywid/2 & y < ywid/2;
%
% % options
% x = CheckOption('x', [], varargin{:}); % default = -N/2:N/2-1 (origin = center)
% y = CheckOption('y', [], varargin{:}); % default = -N/2:N/2-1 (origin = center)
%
% return:
%    Imcrop is the cropped Im
%    bMaskcrop is the cropped bMask
%
% see also PadImArray

% options
x = CheckOption('x', [], varargin{:}); % default = -N/2:N/2-1 (origin = center)
y = CheckOption('y', [], varargin{:}); % default = -N/2:N/2-1 (origin = center)

% center to the nearest pixel
if isempty(xycent), xycent = [0 0]; end

% use CreateGrid instead for consistency with other code
% x = round(nimx/2+1 + xycent(1)) + (-ww/2:ww/2-1)';
% y = round(nimy/2+1 + xycent(2)) + (-ww/2:ww/2-1)';
if isempty(x) || isempty(y),
    [x, y] = CreateGrid(Im); % origin at center
elseif isscalar(x)
    [nr, nc] = size(Im);
    x = x + (1:nc)' - 1;
    y = y + (1:nr)' - 1;
end

x = x - xycent(1);
y = y - xycent(2);

% % make sure cropping does not go outside image boundary
% x = x( x >= 1 & x <= nimx );
% y = y( y >= 1 & y <= nimy );

%
ixuse = x >= -xwid/2 & x < xwid/2;
iyuse = y >= -ywid/2 & y < ywid/2;
Imcrop = Im(iyuse, ixuse);

if all(size(bMask) == size(Im)),
    bMaskcrop = bMask(iyuse, ixuse);
else
    bMaskcrop = [];
end

end % CropImage