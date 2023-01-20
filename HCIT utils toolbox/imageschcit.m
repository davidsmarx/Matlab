
function hout = imageschcit(varargin)
% hh = imageschcit(Im)
% hh = imageschcit(x, y, Im)
% hh = imageschcit(x, y, Im, axes property, value, ...)
%
% Im can be an array or a fits filename
% x, y can be vectors size(Im), or offset scalars
%
% rgb images are size (:,:,3)
% 
% hout return value is handle to image object
%
% other options:
%    'hax', hax % axes handle to put image (default = new axes)
%    'scale', 'linear' or 'log', % default = 'linear'

cProperties = {};
hax = [];
   
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

    % if hax option is specified, get it and remove from the list
    [hax, cProperties] = CheckOptionPop('hax', [], cProperties{:});
    
    % option 'scale'
    [scale, cProperties] = CheckOptionPop('scale', 'linear', cProperties{:});

else,
    error('usage: hh = imageschcit(x, y, Im);');
end

% if Im is a filename, read the image from file
switch class(Im),
    case 'char'
        [pn, fn, fext] = fileparts(Im);
        switch fext,
            case '.fits',
                % check if primary hdu is empty, and if so, use first image hdu
                finfo = fitsinfo(Im);
                if isempty(finfo.PrimaryData.Size) || isequal(finfo.PrimaryData.Size, [1 1]),
                    Im = fitsread(Im,'image');
                else
                    Im = fitsread(Im);
                end
            otherwise,
                error(['image file type ' fext ' not yet implemented']);
        end % switch file extension
    
    case 'double'
        % do nothing
        
    case 'py.numpy.ndarray'
        Im = double(Im);
        
    otherwise
        % if this fails, then type error
        Im = double(Im);
        
end % switch class(Im)


% create grid if necessary
if isempty(x) || isempty(y),
    [x, y] = CreateGrid(Im);
elseif isscalar(x) || isscalar(y),
    % make zero-offset axis to match python arrays
    sizeIm = size(Im); % could be [nr nc] or [nr nc 3] for rgb
    nr = sizeIm(1);
    nc = sizeIm(2);
    if isscalar(x), x = x + (1:nc)'-1; end
    if isscalar(y), y = y + (1:nr)'-1; end
end

maxIm = max(abs(imag(Im(:))));
maxRe = max(abs(real(Im(:))));
if maxIm > eps*maxRe && maxRe > eps*maxIm,
    [hax, hh] = ImageReIm(x, y, Im);
    
else
    if maxIm <= eps*maxRe, Im = real(Im); end
    if maxRe <= eps*maxIm, Im = imag(Im); end
    
    if ~isempty(hax),
        hh = imagesc(hax, x, y, Im);
    else
        hh = imagesc(x, y, Im);
    end
    
    axis image
    set(gca,'ydir','normal')
    colormap(jet)
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
    colormap(jet)
    axis image
    title('Real')
    
    hax(2) = subplot(1,2,2);
    hh(2) = imagesc(x, y, imag(Im));
    colormap(jet)
    axis image
    title('Imag')
    
    set(hax,'ydir','normal');

    clim1 = get(hax(1),'clim');
    clim2 = get(hax(2),'clim');
    set(hax,'clim',[ min([clim1 clim2]) max([clim1 clim2])])
    
end % ImageReIm