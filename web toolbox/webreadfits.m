function [img, finfo] = webreadfits(webaddr, varargin)
% [img, finfo] = webreadfits(webaddr, varargin)

bDisplay = CheckOption('display', true, varargin{:});

%
options_info = weboptions('ContentType', 'image', 'ContentReader', @fitsinfo);
options_read = weboptions('ContentType', 'image', 'ContentReader', @fitsread);

%
finfo = webread(webaddr, options_info);
img = webread(webaddr, options_read);

%
if bDisplay,
    figure, imageschcit(img), colormap gray
end
