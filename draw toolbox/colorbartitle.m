function htout = colorbartitle(varargin)
% ht = colorbartitle(title)
% ht = colorbartitle(hc, title)
% ht = colorbartitle(...,'property name','property value',...)
%
% title is a string
% hc is a handle to a colorbar, if no handle is given then hc = colorbar;
% creates a colorbar in the current figure
%
% return ht = handle to the title in the colorbar axes

if ishandle(varargin{1}),
    hc = varargin{1};
    titlestring = varargin{2};
    options = {varargin{3:end}};
elseif ischar(varargin{1}),
    hc = colorbar;
    titlestring = varargin{1};
    options = {varargin{2:end}};
else
    error('first argument is not a handle or a string');
end

set(get(hc,'label'),'string',titlestring,'Fontsize',12,options{:})

if nargout > 0, htout = hc; end

