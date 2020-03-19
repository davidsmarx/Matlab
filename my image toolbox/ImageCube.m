function [hfig, hax, sUserData] = ImageCube(imgCube, listI, varargin)
% [hfig, hax, sUserData] = ImageCube(imgCube, listI, varargin)
%
% imgCube (Nslices, nr, nc)
% listI (vector of length Nslices) numeric label for each slice
%
% keyboard commands:
%   'f' next slice
%   'b' previous slice
%   '1' first slice
%   'e' last slice
%
% options:
% fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
% xplot = CheckOption('x', [], varargin{:});
% yplot = CheckOption('y', [], varargin{:});
% xylim = CheckOption('xylim', [], varargin{:});
% clim = CheckOption('clim', [], varargin{:});
% cmap = CheckOption('colormap', 'gray', varargin{:});
% fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});
% hfig = CheckOption('hfig', [], varargin{:});
% hax = CheckOption('hax', [], varargin{:});

% validate inputs
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
xylim = CheckOption('xylim', [], varargin{:});
clim = CheckOption('clim', [], varargin{:});
cmap = CheckOption('colormap', 'gray', varargin{:});
fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});


% validate axis scale
if isempty(xplot) || isempty(yplot),
    [xplot, yplot] = CreateGrid([nr nc]);
end

% create the figure
islinit = 1; % first slice to show
if ~isempty(hfig), figure(hfig); else, hfig = gcf; end 
if ~isempty(hax), axes(hax); end
himage = fImageDisplay(xplot, yplot, squeeze(imgCube(islinit,:,:)));
if isempty(hax), hax = gca; end
if ~isempty(clim), set(gca,'clim',clim), end
if ~isempty(xylim), set(gca,'xlim', xylim, 'ylim', xylim), end
colormap(cmap);
htitle = title(fTitleStr(islinit));

% store the image cube, etc., with the figure object
sUserData = struct(...
    'imgCube', imgCube ...
    ,'Nsl', Nsl ...
    ,'hfig', hfig ...
    ,'hax', hax ...
    ,'himage', himage ...
    ,'htitle', htitle ...
    ,'fTitleStr', fTitleStr ...
    ,'isl', islinit ... % current slice in figure
    ,'xplot', xplot ...
    ,'yplot', yplot ...
    );

set(hfig, 'KeyPressFcn', @KeyPressCallback);
set(hfig, 'UserData', sUserData);

end % main

function KeyPressCallback(hSrc, event)

S = get(hSrc,'UserData');

%disp(event.Key)
switch event.Key
    case 'f'
        %disp('forward');
        %S.isl = min(S.Nsl, S.isl + 1);
        S.isl = mod(S.isl, S.Nsl) + 1;
        
    case 'b'
        %disp('backward');
        %S.isl = max(1, S.isl - 1);
        S.isl = S.isl - 1;
        if S.isl == 0, S.isl = S.Nsl; end
        
    case '1'
        % go to first slice
        S.isl = 1;
        
    case 'e'
        % go to last slice
        S.isl = S.Nsl;
        
    case 'm'
        % save as movie?
        [fn, pn] = uiputfile({'*.avi'});
        a = inputdlg('enter frames per sec');
        fps = str2double(a{1});
        if ~isequal(fn,0),
            v = VideoWriter([pn fn]);
            v.FrameRate = fps; % frames/sec
            open(v);            
        end
        
        % make a movie
        loops = S.Nsl;
        F(loops) = struct('cdata',[],'colormap',[]);
        for isl = 1:loops,
            Img = squeeze(S.imgCube(isl,:,:));
            S.himage.CData = Img;
            S.htitle.String = S.fTitleStr(isl);
            drawnow;
            F(isl) = getframe(S.hfig);
            
            if ~isequal(fn,0),
                writeVideo(v, F(isl));
            end
            
        end
   
        if ~isequal(fn,0),
            close(v);
        end
        
        % playback and save as movie
        hfigm = figure; movie(hfigm, F, 2);
        
        
        
    otherwise
        
end

Img = squeeze(S.imgCube(S.isl,:,:));
S.himage.CData = Img;
S.htitle.String = S.fTitleStr(S.isl);

drawnow;

set(hSrc, 'UserData', S);

end % KeyPressCallback

