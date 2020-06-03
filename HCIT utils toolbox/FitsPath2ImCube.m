function varargout = FitsPath2ImCube(pn, varargin)
% [ImCube, sParms] = FitsPath2ImCube(pn, varargin)
% ImCube = FitsPath2ImCube(pn, options)
%
% simple routine to collect all the fits files in a subdir
% for example, making an image cube out of all the fits files in a
% phase retrieval subdir
%
% options:
% plottype = CheckOption('plottype', 'cube', varargin{:}); % 'spread', 'cube', 'none'
% plotx = CheckOption('x', 0, varargin{:});
% ploty = CheckOption('y', 0, varargin{:});
% xlim = CheckOption('xlim', [], varargin{:});
% ylim = CheckOption('ylim', [], varargin{:});
% hdrkwd = CheckOption('hdrkwd', {'camz'}, varargin{:});
% hdrkwdvalfmt = CheckOption('hdrkwdvalfmt', '%.1f', varargin{:});
% scale = CheckOption('scale', 'linear', varargin{:}); or 'log'
% refImg = CheckOption('refimg', [], varargin{:}); ImCube = ImCube - refImg
% comTitlestr = CheckOption('comtitlestr', '', varargin{:}); % common title to start each titlestr
% 
% output:
% ImCube = [Nimages nr nc]
% sParms = struct(

plottype = CheckOption('plottype', 'cube', varargin{:}); % 'spread'
plotx = CheckOption('x', 0, varargin{:});
ploty = CheckOption('y', 0, varargin{:});
xlim = CheckOption('xlim', [], varargin{:});
ylim = CheckOption('ylim', [], varargin{:});
scale = CheckOption('scale', 'linear', varargin{:});
hdrkwd = CheckOption('hdrkwd', {'camz'}, varargin{:});
hdrkwdvalfmt = CheckOption('hdrkwdvalfmt', '%.1f ', varargin{:});
refImg = CheckOption('refimg', [], varargin{:}); % ImCube = ImCube - refImg
comTitlestr = CheckOption('comtitlestr', '', varargin{:}); % common title to start each titlestr

%

% initialize return vals
hfig = [];
hax = [];
 
listfn = dir(PathTranslator([pn '/*.fits']));
Nf = length(listfn);

% read one image to get image size and extension
finfo = fitsinfo(PathTranslator([pn '/' listfn(1).name]));
if isempty(finfo.PrimaryData.Size),
    sExtension = 'image';
else
    sExtension = 'primary';
end
imtmp = fitsread(PathTranslator([pn '/' listfn(1).name]),sExtension);

[Nr, Nc] = size(imtmp);

% check header keyword, make sure it's cell
if ischar(hdrkwd), hdrkwd = {hdrkwd}; end

% allocate
ImCube = zeros([Nf Nr Nc]);
hdrkwdval = zeros(Nf,length(hdrkwd));

for ii = 1:Nf

    imtmp = fitsread(PathTranslator([pn '/' listfn(ii).name]),sExtension);
    finfo = fitsinfo(PathTranslator([pn '/' listfn(ii).name]));

    if ~isempty(refImg),
        imtmp = imtmp - refImg;
    end
    
    ImCube(ii,:,:) = imtmp;

    % this way, FitsGetKeywordVal always returns a cell array
    % even if only one hdrkwd
    ctmp = FitsGetKeywordVal(finfo.PrimaryData.Keywords,hdrkwd);
    hdrkwdval(ii,:) = [ctmp{:}];
    
end

% apply scale (stretch)
switch lower(scale),
    case 'linear',
        % do nothing
    case 'log'
        ImCube = log10(ImCube);
        
    otherwise,
        error(['unknown scale: ' scale]);
end


switch lower(plottype),
    case 'cube',
        if ~isempty(hdrkwd)
            fTitleStr = @(isl) [[comTitlestr '#' num2str(isl)] join(string(hdrkwd), ', ') sprintf(hdrkwdvalfmt, hdrkwdval(isl,:))];
        else
            fTitleStr = @(isl) ['# ' num2str(isl)];
            hdrkwdval = 1:Nf; % just to have labels for the image cube slices
        end
        
        figure,        
        [hfig, hax, sUserData] = ImageCube(ImCube, hdrkwdval, ...
            'fTitleStr', fTitleStr, ...
            'x', plotx, 'y', ploty);
        
    case 'spread'
        figure
        [hfig, hax] = PlotSpread;

    case 'none',
        % no dispaly, do nothing
        
    otherwise,
        error(['unknown plottype ' plottype]);
end

% return values, depends on options
if nargout == 0,
    varargout = {};
else
    % 
    sParms = struct(...
        'hfig', hfig ...
        ,'hax', hax ...
        ,'listfn', listfn ...
        );
    for ikwd = 1:length(hdrkwd) % hdrkwd is a cell array
        sParms.(hdrkwd{ikwd}) = hdrkwdval(:,ikwd);
    end

    varargout = {ImCube, sParms};
end


return % end of main

%%%%%%%%%%%%%% nested plot functions
    % nested function to create panoramic spread plot has all the options
    % already in scope
    function [hfig, hax] = PlotSpread
        hfig = gcf;
        
        % x, y axes
        if isscalar(plotx),
            x = (plotx:(Nc-plotx-1))';
        else
            x = plotx;
        end
        if isscalar(ploty),
            y = (ploty:(Nr-ploty-1))';
        else
            y = ploty;
        end
        
        % apply xlim, ylim directly to ImCube
        if ~isempty(xlim)
            ix = x>=xlim(1) & x<=xlim(2);
            ImCube = ImCube(:,:,ix);
            x = x(ix);
        end
        if ~isempty(ylim)
            iy = (y>=ylim(1) & y<=ylim(2));
            ImCube = ImCube(:,iy,:);
            y = y(iy);
        end
        
        % unfold cube
        nr = length(y); nc = length(x);
        img = zeros(nr, Nf*nc);
        for isl = 1:Nf,
            img(:, (isl-1)*nc + (1:nc)) = squeeze(ImCube(isl,:,:));
        end %
        imageschcit(img)
        hax = gca;
        title(pwd2titlestr(pn),'fontsize',14)
        
    end % PlotSpread

end % main

