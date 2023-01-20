function [hfig_out, hax_out, sUserData_out] = ImageCube(imgCube, listI, varargin)
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
%   'm' make a movie stepping through the slices and save as .avi
%
% options:
% fImageDisplay = CheckOption('fImageDisplay', @imageschcit, varargin{:});
% xplot = CheckOption('x', [], varargin{:});
% yplot = CheckOption('y', [], varargin{:});
% xylim = CheckOption('xylim', [], varargin{:});
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

imcube_size = size(imgCube); % allow for rgb
Nsl = imcube_size(1);
nr = imcube_size(2);
nc = imcube_size(3);

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
    ,'isl', islinit ... % current slice in figure
    ,'xplot', xplot ...
    ,'yplot', yplot ...
    ,'multi_key_seq', [] ...
    );
sUserData.fTitleStr = fTitleStr; % allows for fTitleStr is a cell array

set(hfig, 'KeyPressFcn', @KeyPressCallback);
set(hfig, 'UserData', sUserData);

if nargout > 0
    [hfig_out, hax_out, sUserData_out] = deal(hfig, hax, sUserData);
end

end % main

function KeyPressCallback(hSrc, event)

S = get(hSrc,'UserData');

% are we in the middle of a multi-key sequence?
if ~isempty(S.multi_key_seq),
    S.multi_key_seq = [S.multi_key_seq event.Key];    
end

%disp(event.Key)
% check for 'shift', 'control', etc
if any(strcmp(event.Modifier,'shift'))
    switch event.Key
        case 'g'
            % export to GIF, requires ver 2022a ?
            [fn, pn] = uiputfile({'*.gif'});
            if isequal(fn, 0)
                % user pressed Cancel
                return
            end
            
            % loop through frames
            loops = S.Nsl;
            for isl = 1:loops,
                Img = squeeze(S.imgCube(isl,:,:));
                S.himage.CData = Img;
                S.htitle.String = S.fTitleStr(isl);
                drawnow;
                exportgraphics(gcf, fullfile(pn, fn), 'Append', true);
            end
        otherwise
            
    end % switch 'shift' event.Key
    
%elseif isequal(event.Modifier,'control')
    
    
else
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
            if isempty(S.multi_key_seq),
                % go to first slice
                S.isl = 1;
            end
            
        case 'e'
            % go to last slice
            S.isl = S.Nsl;
            
        case 'g'
            % goto frame number
            % command is g number g
            % start or end multi-key sequence
            if isempty(S.multi_key_seq),
                % start new multi-key sequence
                S.multi_key_seq = 'g';
            else
                % check that multi key sequence starts with 'g'
                if S.multi_key_seq(1) == 'g',
                    % assumes all keys in the middle are numbers
                    isl = str2double(S.multi_key_seq(2:end-1));
                    if ~isnan(isl),
                        S.isl = isl;
                    end
                end
                % reset
                S.multi_key_seq = '';
                
            end
            
        case 'm'
            % save as movie?
            [fn, pn] = uiputfile({'*.avi'});
            if isequal(fn, 0)
                % user pressed Cancel
                return
            end
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
            
    end % switch event.Key
    
end % if 'shift', etc

if S.isl >=1 && S.isl <= S.Nsl,
    Img = squeeze(S.imgCube(S.isl,:,:));
    S.himage.CData = Img;
    S.htitle.String = S.fTitleStr(S.isl);
end

drawnow;

set(hSrc, 'UserData', S);

end % KeyPressCallback

