function [img, finfo] = webreadfits(webaddr, varargin)
% [img, finfo] = webreadfits(webaddr, varargin)
%
% CheckOption('display', true, varargin{:});
% CheckOption('displayscale', 'log', varargin{:});

bDisplay = CheckOption('display', true, varargin{:});
scaleDisplay = CheckOption('displayscale', 'log', varargin{:});
extname = CheckOption('extname', 'image', varargin{:}); % 'primary', 'image', 

%
options_info = weboptions('ContentType', 'image', 'ContentReader', @fitsinfo);
%options_read = weboptions('ContentType', 'image', 'ContentReader', @fitsread);
options_read = weboptions('ContentType', 'image', 'ContentReader', @ReadFitsImage);

%
finfo = webread(webaddr, options_info);
img = webread(webaddr, options_read);

%
if bDisplay,
    switch scaleDisplay,
        case 'log'
            fScale = @(img) logImage(img);
        otherwise
            fScale = @(img) img;
    end

    figure, imageschcit(fScale(img)), colormap gray
end

    function img = ReadFitsImage(fn)
        
        img = fitsread(fn, extname);
        
    end

end