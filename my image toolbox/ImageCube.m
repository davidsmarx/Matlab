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
% clim = CheckOption('clim', [], varargin{:});
% cmap = CheckOption('colormap', 'gray', varargin{:});
% fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});
% hfig = CheckOption('hfig', [], varargin{:});
% hax = CheckOption('hax', [], varargin{:});

% check options
hfig = CheckOption('hfig', [], varargin{:});
hax = CheckOption('hax', [], varargin{:});
fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
clim = CheckOption('clim', [], varargin{:});
xylim = CheckOption('xylim', [], varargin{:});
cmap = CheckOption('colormap', 'gray', varargin{:});
fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});

% validate inputs
[Nsl, nr, nc] = size(imgCube);
if length(listI) ~= Nsl,
    error('number of image slices not consistent');
end

% create the figure
islinit = 1; % first slice to show
if ~isempty(hfig), figure(hfig); else, hfig = gcf; end 
if ~isempty(hax), axes(hax); end
himage = fImageDisplay(squeeze(imgCube(islinit,:,:)));
if isempty(hax), hax = gca; end
if ~isempty(clim), set(gca,'clim',clim), end
if ~isempty(xylim), set(gca,'xlim', xylim, 'ylim', xylim), end
colormap(cmap);
htitle = title(fTitleStr(islinit));

% store the image cube, etc., with the figure object
sUserData = struct(...
    'imgCube', imgCube ...
    ,'Nsl', Nsl ...
    ,'listSlLabel', listI ... % not used
    ,'hax', hax ...
    ,'himage', himage ...
    ,'htitle', htitle ...
    ,'fTitleStr', fTitleStr ...
    ,'isl', islinit ... % current slice in figure
    );

set(hfig, 'KeyPressFcn', @KeyPressCallback);
set(hfig, 'UserData', sUserData);


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
                
            otherwise
                
        end 
        
        Img = squeeze(S.imgCube(S.isl,:,:));
        S.himage.CData = Img;
        S.htitle.String = S.fTitleStr(S.isl);
        
        drawnow;
        
        set(hSrc, 'UserData', S);
        
    end % KeyPressCallback

end % main

