function [hfig, sUserData] = ImageCube(imgCube, listI, varargin)
% [hfig, sUserData] = ImageCube(imgCube, listI, varargin)
%
% imgCube (Nslices, nr, nc)
% listI (vector of length Nslices) numeric label for each slice
%
% options:
% fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
% clim = CheckOption('clim', [], varargin{:});
% cmap = CheckOption('colormap', 'gray', varargin{:});
% fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});

% check options
hfig = CheckOption('hfig', [], varargin{:});
fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
clim = CheckOption('clim', [], varargin{:});
cmap = CheckOption('colormap', 'gray', varargin{:});
fTitleStr = CheckOption('fTitleStr', @(isl) ['slice #' num2str(isl) '; Label ' num2str(listI(isl))], varargin{:});

% validate inputs
[Nsl, nr, nc] = size(imgCube);
if length(listI) ~= Nsl,
    error('number of image slices not consistent');
end

% create the figure
islinit = 1; % first slice to show
if isempty(hfig), hfig = figure; else, figure(hfig), end
himage = fImageDisplay(squeeze(imgCube(islinit,:,:)));
if ~isempty(clim), set(gca,'clim',clim), end
colormap(cmap);
htitle = title(fTitleStr(islinit));

% store the image cube, etc., with the figure object
sUserData = struct(...
    'imgCube', imgCube ...
    ,'Nsl', Nsl ...
    ,'listSlLabel', listI ... % not used
    ,'hax', gca ...
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
                S.isl = min(S.Nsl, S.isl + 1);
                
            case 'b'
                %disp('backward');
                S.isl = max(1, S.isl - 1);

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

