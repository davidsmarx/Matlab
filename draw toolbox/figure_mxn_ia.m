function hf = figure_mxn_ia(m,n,varargin)
% hfig = figure_mxn(m,n,varargin)
% 
% for creating m rows by n columns of subplots

hfig = figure;
newpos = get(hfig,'position').*[1 1 n m];

% check against screen 
screenlims = get(0,'ScreenSize');
if newpos(4) > screenlims(4),
    newpos(3:4) = (screenlims(4)./newpos(4)).*newpos(3:4);
end
if newpos(2)+newpos(4) > screenlims(4)-80,
    newpos(2) = screenlims(4) - newpos(4) - 80;
end

set(hfig,'position',newpos)

if nargin > 2,
    set(hfig,varargin{:})
end

if nargout >= 1,
    hf = hfig;
end