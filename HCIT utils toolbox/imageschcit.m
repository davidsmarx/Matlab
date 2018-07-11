function hout = imageschcit(varargin)
% hh = imageschcit(Im)
% hh = imageschcit(x, y, Im)
% hh = imageschcit(x, y, Im, axes property, value, ...)

cProperties = {};
        
if nargin == 1,
    Im = varargin{1};
    nn = size(Im); % in case Im is more than 2-d, e.g. rgb
    nr = nn(1);
    nc = nn(2);
    x = (1:nc)';
    y = (1:nr)';
    
elseif nargin == 3,
    [x, y, Im] = deal(varargin{:});
    
elseif nargin > 3,
    [x, y, Im] = deal(varargin{1:3});
    cProperties = {varargin{4:end}};
    
else,
    error('usage: hh = imageschcit(x, y, Im);');
end

hh = imagesc(x, y, Im);
axis image
set(gca,'ydir','normal')

if ~isempty(cProperties),
    set(gca, cProperties{:});
end


% return value?
if nargout > 0,
    hout = hh;
end