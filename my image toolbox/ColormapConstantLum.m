function cmapout = ColormapConstantLum(cmap)
% cmapout = ColormapConstantLum
% cmapout = ColormapConstantLum(cmap)
%
% convert cmap to a colormap with constant luminance
% default cmap = colormap(cool)

if nargin < 1 || isempty(cmap),
    cmap = colormap(cool);
end

mycmap = cmap ./ rgb2gray(cmap);
cmapout = mycmap ./ max(mycmap(:));

return

