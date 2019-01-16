function hout = imageschcit(varargin)
% hh = imageschcit(Im)
% hh = imageschcit(x, y, Im)
% hh = imageschcit(x, y, Im, axes property, value, ...)
%
% Im can be an array or a fits filename
% x, y can be vectors size(Im), or offset scalars

cProperties = {};
   
x = []; y = [];
if nargin == 1,
    Im = varargin{1};
    %     nn = size(Im); % in case Im is more than 2-d, e.g. rgb
    %     nr = nn(1);
    %     nc = nn(2);
    %     x = (1:nc)';
    %     y = (1:nr)';
    
elseif nargin == 3,
    [x, y, Im] = deal(varargin{:});
    
elseif nargin > 3,
    [x, y, Im] = deal(varargin{1:3});
    cProperties = {varargin{4:end}};
    
else,
    error('usage: hh = imageschcit(x, y, Im);');
end

% if Im is a filename, read the image from file
if ischar(Im),
    [pn, fn, fext] = fileparts(Im);
    switch fext,
        case '.fits',
            % check if primary hdu is empty, and if so, use first image hdu
            finfo = fitsinfo(Im);
            if isempty(finfo.PrimaryData.Size),
                Im = fitsread(Im,'image');
            else
                Im = fitsread(Im);
            end
        otherwise,
            error(['image file type ' fext ' not yet implemented']);
    end
end

% create grid if necessary
if isempty(x) || isempty(y),
    [x, y] = CreateGrid(Im);
elseif isscalar(x) || isscalar(y),
    % make zero-offset axis to match python arrays
    [nr, nc] = size(Im);
    if isscalar(x), x = x + (1:nc)'-1; end
    if isscalar(y), y = y + (1:nr)'-1; end
end

maxIm = max(abs(imag(Im(:))));
maxRe = max(abs(real(Im(:))));
if maxIm > eps*maxRe && maxRe > eps*maxIm,
    [hax, hh] = ImageReIm(x, y, Im);
    
else
    hh = imagesc(x, y, Im);
    axis image
    set(gca,'ydir','normal')
    hax = gca;
end

if ~isempty(cProperties),
    set(hax, cProperties{:});
end


% return value?
if nargout > 0,
    hout = hh;
end

end % main

function [hax, hh] = ImageReIm(x, y, Im)

    hax(1) = subplot(1,2,1);
    hh(1) = imagesc(x, y, real(Im));
    axis image
    
    hax(2) = subplot(1,2,2);
    hh(2) = imagesc(x, y, imag(Im));
    axis image
    
    set(hax,'ydir','normal');

    clim1 = get(hax(1),'clim');
    clim2 = get(hax(2),'clim');
    set(hax,'clim',[ min([clim1 clim2]) max([clim1 clim2])])
    
end % ImageReIm