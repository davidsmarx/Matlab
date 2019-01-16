function hf = figure_mxn(varargin)
% hfig = figure_mxn(m,n,listProperties)
% hfig = figure_mxn([m, n], listProperties)
% figure_mxn(hfig, ...)
% 
% for creating m rows by n columns of subplots

if isa(varargin{1},'matlab.ui.Figure'),
    hfig = varargin{1};
    m = varargin{2};
    n = varargin{3};
    listProps = varargin(4:end);
elseif isscalar(varargin{1}),
    m = varargin{1};
    n = varargin{2};
    listProps = varargin(3:end);
    hfig = figure;
else
    m = varargin{1}(1);
    n = varargin{1}(2);
    listProps = varargin(2:end);
    hfig = figure;
end

% standard figure size
dx = 560; dy = 420;
newpos = get(hfig,'position').*[1 1 0 0] + [0 0 dx*n dy*m];

% check against screen 
screenlims = get(0,'ScreenSize');
if newpos(4) > screenlims(4),
    newpos(3:4) = (screenlims(4)./newpos(4)).*newpos(3:4);
end
if newpos(2)+newpos(4) > screenlims(4)-80,
    newpos(2) = screenlims(4) - newpos(4) - 80;
end

set(hfig,'position',newpos);

if ~isempty(listProps),
    set(hfig,listProps{:})
end

if nargout >= 1,
    hf = hfig;
end