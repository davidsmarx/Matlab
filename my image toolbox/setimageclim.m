function setimageclim(h,clim)
% setimageclim(h,clim)
%
% h = handle to a figure or axes with an image
% clim = new limits, or 'auto' to auto set (default = auto)
%   for the image display

if nargin == 0, h = gcf; end

switch get(h,'Type')
    case 'figure'
        h_fig = h;
        figure(h_fig);
        h_iaxis = gca;
    case 'axes'
        h_fig = get(h,'Parent');
        h_iaxis = h;
        figure(h_fig);
    otherwise,
        error(['handle type ' get(h,'Type') ' not supported']);
end

% get the handle to the image
hc = get(h_iaxis,'children');
h_image = hc( strmatch('image',get(hc,'type')) );
if isempty(h_image), error('no image found'); end

% check for auto option
if ~exist('clim','var'), clim = 'auto'; end
if strcmp(clim,'auto')
    clim = autoclim(h_image); % returns empty matrix if fails
end

if ~isempty(clim),
    set(h_iaxis,'clim',clim);
    set(h_image(isprop(h_image,'cdatamapping')),'cdatamapping','scaled');
else,
    warning('could not reset clim');
end

return

function clim = autoclim(h_image)

cdata = get(h_image,'cdata');
[chist, x] = hist(double(cdata(:)),256);

cumdist = cumsum(chist)./sum(chist);

xsmall = [x(1) x(cumdist < 0.05)];
% check that not all the pixels are zero (black background with a narrow
% peak)
if cumdist(1) > 0.99,
    warning('histogram too narrow, could not auto clim');
    clim = [];
else, % normal auto clim
    xbig = [x(cumdist > 0.99) x(end)];
    clim = [xsmall(end) xbig(1)];
end

return

