function h = imageopd(x,y,opd,varargin)
% h = imageopd(x,y,opd,options)
%
% options:
%    'title'
%    '
flipmap = @(a) (fliplr(rot90(a)));

hh = imagesc(y,x,flipmap(opd));
%hh = pcolor(y,x,flipmap(opd)); shading interp;

% if any(isnan(opd(:))),
%     % adjust colormap
%     cmap = colormap;
%     colormap([1 1 1; cmap]);
% end

% get axis handle
ha = get(hh,'Parent');
set(ha,'YDir','normal');
axis image
xlabel('Y (deg)'), ylabel('X (deg)')
hc = colorbar;

% check for options
for ii = 1:length(varargin),
    switch lower(varargin{ii}),
        case 'title'
            title(varargin{ii+1});
            ii = ii+1;
        case 'colorbartitle'
            set(get(hc,'Title'),'String',varargin{ii+1});
            ii = ii+1;
        otherwise,
            % nothing
    end
end

if nargout > 0,
    h = hh;
end