function [hfig_out, hax_out, sUserData_out] = ImageCubeSpread(imgCube, listI, varargin)
% [hfig, hax, sUserData] = ImageCubeSpread(imgCube, listI, varargin)
%
% imgCube (Nslices, nr, nc)
% listI (vector of length Nslices) numeric label for each slice
%
% keyboard commands:
%   'f' next slice
%   'b' previous slice
%   '1' first slice
%   'e' last slice
%   'm' make a movie stepping through the slices and save as .avi
%
% options:
% fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
% xplot = CheckOption('x', [], varargin{:});
% yplot = CheckOption('y', [], varargin{:});
% xlim = CheckOption('xlim', [], varargin{:}); % xlim for each slice
% ylim = CheckOption('ylim', [], varargin{:}); % ylim for each slice
% clim = CheckOption('clim', [], varargin{:});
% cmap = CheckOption('colormap', 'gray', varargin{:});
% fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});
% hfig = CheckOption('hfig', [], varargin{:}); % put image in hfig, or gcf
% hax = CheckOption('hax', [], varargin{:});
%
% fTitleStr can also be a cell array of strings

% validate inputs
if ischar(imgCube)
    imgCube = shiftdim(fitsread(imgCube),2);
end

[Nsl, nr, nc] = size(imgCube);
if nargin == 1,
    listI = 1:Nsl;
end
if length(listI) ~= Nsl,
    error('number of image slices not consistent');
end

% check options
hfig = CheckOption('hfig', [], varargin{:});
hax = CheckOption('hax', [], varargin{:});
fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
xplot = CheckOption('x', [], varargin{:});
yplot = CheckOption('y', [], varargin{:});
xlim = CheckOption('xlim', [], varargin{:});
ylim = CheckOption('ylim', [], varargin{:});
clim = CheckOption('clim', [], varargin{:});
cmap = CheckOption('colormap', 'gray', varargin{:});
fTitleStr = CheckOption('fTitleStr', @(isl) ['#' num2str(listI(isl))], varargin{:});


% validate axis scale
if isempty(xplot) || isempty(yplot),
    [xplot, yplot] = CreateGrid([nr nc]);
end

% create the figure
if ~isempty(hfig), figure(hfig); else, hfig = gcf; end 
%if ~isempty(hax), axes(hax); end
%himage = fImageDisplay(xplot, yplot, squeeze(imgCube(islinit,:,:)));
% if isempty(hax), hax = gca; end
% if ~isempty(clim), set(gca,'clim',clim), end
% if ~isempty(xylim), set(gca,'xlim', xylim, 'ylim', xylim), end
% colormap(cmap);
% htitle = title(fTitleStr(islinit));

% % store the image cube, etc., with the figure object
% sUserData = struct(...
%     'imgCube', imgCube ...
%     ,'Nsl', Nsl ...
%     ,'hfig', hfig ...
%     ,'hax', hax ...
%     ,'himage', himage ...
%     ,'htitle', htitle ...
%     ,'isl', islinit ... % current slice in figure
%     ,'xplot', xplot ...
%     ,'yplot', yplot ...
%     ,'multi_key_seq', [] ...
%     );
% sUserData.fTitleStr = fTitleStr; % allows for fTitleStr is a cell array

% x, y axes
if isscalar(xplot),
    x = (xplot:(nc-xplot-1))';
else
    x = xplot;
end
if isscalar(yplot),
    y = (yplot:(nr-yplot-1))';
else
    y = yplot;
end

% apply xlim, ylim directly to imgCube
if ~isempty(xlim)
    ix = x>=xlim(1) & x<=xlim(2);
    imgCube = imgCube(:,:,ix);
    x = x(ix);
end
if ~isempty(ylim)
    iy = (y>=ylim(1) & y<=ylim(2));
    imgCube = imgCube(:,iy,:);
    y = y(iy);
end

% unfold cube
nr = length(y); nc = length(x);
img = zeros(nr, Nsl*nc);
for isl = 1:Nsl,
    img(:, (isl-1)*nc + (1:nc)) = squeeze(imgCube(isl,:,:));
end %
himg = imageschcit(1:Nsl*nc, 1:nr, img);
hax = gca;
hax.XTick = [[1:nc:Nsl*nc] Nsl*nc];
hax.XTickLabel = cat(2, repmat({'1'}, [1 Nsl]), {num2str(nc)});

% titles
htitle(1) = title(fTitleStr(1));
posTitle = htitle(1).Position;
htitle(1).Position = [0.5*nc posTitle(2:3)];
for isl = 2:Nsl,
    htitle(isl) = copy(htitle(1));
    htitle(isl).Parent = hax;
    htitle(isl).String = fTitleStr(isl);
    htitle(isl).Position = [(isl-1)*nc+0.5*nc posTitle(2:3)];
end

% resize figure
hfig.Position(3) = Nsl*hfig.Position(4);
screensize = get(0, "ScreenSize");
if hfig.Position(3) > screensize(3)
    a = 0.95*screensize(3)/hfig.Position(3);
    hfig.Position(3:4) = a*hfig.Position(3:4);    
end

% clim and cmap
if ~isempty(clim), set(gca,'clim',clim), end
colormap(cmap);


% return values
if nargout > 0,
    hfig_out = hfig;
    hax_out = hax;
    sUserData_out = struct('himg', himg);
end