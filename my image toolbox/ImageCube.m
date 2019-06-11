function [hfig, sUserData] = ImageCube(imgCube, listI, varargin)

% libFun = struct(...
%     'KeyPressCallback', @KeyPressCallback ...
%     );

fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
clim = CheckOption('clim', [], varargin{:});
cmap = CheckOption('colormap', 'gray', varargin{:});
fTitleStr = CheckOption('fTitleStr', @(ii) ['ii = ' num2str(ii)], varargin{:});

hfig = figure;
himage = fImageDisplay(squeeze(imgCube(1,:,:)));
if ~isempty(clim), set(gca,'clim',clim), end
colormap(cmap);
htitle = title(fTitleStr(1));

sUserData = struct(...
    'imgCube', imgCube ...
    ,'listii', listI ...
    ,'hax', gca ...
    ,'himage', himage ...
    ,'htitle', htitle ...
    ,'fTitleStr', fTitleStr ...
    ,'ii', 1 ...
    );

set(hfig, 'KeyPressFcn', @KeyPressCallback);
set(hfig, 'UserData', sUserData);


    function KeyPressCallback(hSrc, event)
        
        S = get(hSrc,'UserData');
        
        %disp(event.Key)
        switch event.Key
            case 'f'
                %disp('forward');
                S.ii = min(length(S.listii), S.ii + 1);
                
            case 'b'
                %disp('backward');
                S.ii = max(1, S.ii - 1);
                
            otherwise
                
        end 
        
        Img = squeeze(S.imgCube(S.ii,:,:));
        S.himage.CData = Img;
        S.htitle.String = S.fTitleStr(S.ii);
        
        drawnow;
        
        set(hSrc, 'UserData', S);
        
    end % KeyPressCallback

end % main

