function hf = figure_1x2(varargin)
% hfig = figure_1x2(varargin)

hfig = figure;
set(hfig,'position',get(hfig,'position').*[1 1 2 1]);

if nargin > 2,
    set(hfig,varargin{:})
end

if nargout >= 1,
    hf = hfig;
end