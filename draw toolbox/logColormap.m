function cmap_log = logColormap(varargin)
% cmap_log = logColormap(cmap, varargin)
%
% log image ds9 style by mapping the colormap
%
% from http://ds9.si.edu/doc/ref/how.html
%
% CheckOption('alpha', 1000, varargin{:});
% CheckOption('cmap', colormap('gray'), varargin{:}); % can be array or string

% default value
cmap = colormap('gray');

% 
if length(varargin) == 1
    cmap = varargin{1};
end

cmap = CheckOption('cmap', cmap, varargin{:}); % can be array or string
alpha = CheckOption('alpha', 1000, varargin{:});

if ischar(cmap)
    cmap = colormap(cmap);
end

% check cmap is already scaled 0 to 1
if min(cmap(:)) < 0 || max(cmap(:)) > 1,
    error('cmap must be 0 to 1')
end

cmap_log = log10(alpha*cmap + 1)./log10(alpha + 1);

end

