function [Imcrop, bMaskcrop] = CropImage(Im, bMask, xycent, xwid, ywid)
% [Imcrop, bMaskcrop] = CropImage(Im, bMask, xycent, xwid, ywid)
% simple square cropping of the image and its bMask
%
% xycent = [xc yc] in pixels, with [0, 0] at the center of the image array
% xwid, ywid are in pixels
%
% crop to square array of width = max([xwid ywid])
% return:
%    Imcrop is the cropped Im
%    bMaskcrop is the cropped bMask
%
% see also PadImArray

[nimy, nimx] = size(Im);

ww = max([xwid ywid]);
ww = 2*ceil(ww/2); % make sure image size is even

% center to the nearest pixel
if isempty(xycent), xycent = [0 0]; end
x = round(nimx/2+1 + xycent(1)) + (-ww/2:ww/2-1)';
y = round(nimy/2+1 + xycent(2)) + (-ww/2:ww/2-1)';

% make sure cropping does not go outside image boundary
x = x( x >= 1 & x <= nimx );
y = y( y >= 1 & y <= nimy );

%
Imcrop = Im(y,x);

if all(size(bMask) == size(Im)),
    bMaskcrop = bMask(y,x);
else
    bMaskcrop = [];
end

end % CropImage