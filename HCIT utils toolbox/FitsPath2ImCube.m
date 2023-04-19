function varargout = FitsPath2ImCube(pn, varargin)
% [ImCube, sParms] = FitsPath2ImCube(pn, varargin)
% ImCube = FitsPath2ImCube(pn, options)
%
% simple routine to collect all the fits files in a subdir
% for example, making an image cube out of all the fits files in a
% phase retrieval subdir
%
% input:
%   pn
%      string = a valid path, then all pn/*.fits files are used
%      struct = result of dir('...'), then all files in the struct array
%      
%
% options:
% plottype = CheckOption('plottype', 'cube', varargin{:}); % 'spread', 'cube', 'none'
% plotx = CheckOption('x', 0, varargin{:}); % default = 0-offset
% ploty = CheckOption('y', 0, varargin{:}); % default = 0-offset
% xlim = CheckOption('xlim', [], varargin{:});
% ylim = CheckOption('ylim', [], varargin{:});
% hdrkwd = CheckOption('hdrkwd', {'camz'}, varargin{:});
% hdrkwdvalfmt = CheckOption('hdrkwdvalfmt', '%.1f', varargin{:});
% scale = CheckOption('scale', 'linear', varargin{:}); or 'log'
% refImg = CheckOption('refimg', [], varargin{:}); ImCube = ImCube - refImg
% comTitlestr = CheckOption('comtitlestr', '', varargin{:}); % common title to start each titlestr
% clim = CheckOption('clim', [], varargin{:});
% hifg = CheckOption('hfig', [], varargin{:}); % if empty, and plottype ~= 'none', FitsPath2ImCube creates a new figure
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
hdrkwdvalfmt = CheckOption('hdrkwdvalfmt', '%.3f ', varargin{:});
refImg = CheckOption('refimg', [], varargin{:}); % ImCube = ImCube - refImg
comTitlestr = CheckOption('comtitlestr', '', varargin{:}); % common title to start each titlestr
clim = CheckOption('clim', [], varargin{:});
hfig = CheckOption('hfig', [], varargin{:}); % if empty, and plottype ~= 'none', FitsPath2ImCube creates a new figure

%

% initialize return vals
hax = [];
 
switch class(pn)
    case 'char'
        listfn = dir(PathTranslator([pn '/*.fits']));
    case 'struct'
        listfn = pn; % assumes pn is return value from dir()
        pn = listfn(1).folder;
    otherwise
        error(['pn is class ' class(pn)]);
end
        
Nf = length(listfn);

% read one image to get image size and extension
finfo = fitsinfo(PathTranslator([pn '/' listfn(1).name]));
if isempty(finfo.PrimaryData.Size),
    sExtension = 'image';
else
    sExtension = 'primary';
end
%imtmp = fitsread(PathTranslator([pn '/' listfn(1).name]),sExtension);

% check header keyword, make sure it's cell
if ischar(hdrkwd), hdrkwd = {hdrkwd}; end

% allocate
hdrkwdval = zeros(Nf,length(hdrkwd));
[nr_ii, nc_ii] = deal(zeros(Nf,1));
img_ii = cell(Nf,1);
list_ii_skip = [];
for ii = 1:Nf

    img_tmp = fitsread(PathTranslator([pn '/' listfn(ii).name]),sExtension);
    finfo = fitsinfo(PathTranslator([pn '/' listfn(ii).name]));
    
    % can only handle 2-d images
    if ~isequal(ndims(img_tmp), 2)
        list_ii_skip = [list_ii_skip ii];
        continue
    end
    
    %
    img_ii{ii} = img_tmp;
    [nr_ii(ii), nc_ii(ii)] = size(img_ii{ii});
    
    % this way, FitsGetKeywordVal always returns a cell array
    % even if only one hdrkwd
    ctmp = FitsGetKeywordVal(finfo.PrimaryData.Keywords,hdrkwd);
    %     if any(cellfun(@isempty,ctmp)),
    %         ctmp{isempty(ctmp)} = nan;
    %     end
    hdrkwdval(ii,:) = [ctmp{:}];
    %hdrkwdval(ii,:) = 1;
    
end
% remove skipped
if ~isempty(list_ii_skip)
    img_ii(list_ii_skip) = [];
    hdrkwdval(list_ii_skip,:) = [];
    nr_ii(list_ii_skip) = [];
    nc_ii(list_ii_skip) = [];
    listfn(list_ii_skip) = [];
end
Nf = length(img_ii);

% build image cube
Nr = max(nr_ii);
Nc = max(nc_ii);
ImCube = zeros([Nf Nr Nc]);
for ii = 1:Nf
    ImCube(ii,:,:) = PadImArray(img_ii{ii}, [Nr Nc]);
end

% if subtract reference image
if ~isempty(refImg),
    if ~isequal(size(refImg), [Nr Nc])
        refImg = PadImArray(refImg, [Nr Nc]);
    end    
    ImCube = ImCube - shiftdim(repmat(refImg, [1 1 Nf]));
end

% apply scale (stretch)
switch lower(scale),
    case 'linear',
        % do nothing
        ImCube_disp = ImCube;
    case 'log'
        %ImCube = log10(ImCube);
        % scale to 0..1
        ImCube_disp = logImage(ImCube);
        
    otherwise,
        error(['unknown scale: ' scale]);
end

% title strings
if ~isempty(hdrkwd)
    fTitleStr = @(isl) [[comTitlestr '#' num2str(isl)] join(string(hdrkwd), ', ') sprintf(hdrkwdvalfmt, hdrkwdval(isl,:))];
else
    fTitleStr = @(isl) ['# ' num2str(isl)];
    hdrkwdval = 1:Nf; % just to have labels for the image cube slices
end

switch lower(plottype),
    case 'cube',
        
        if isempty(hfig), hfig = figure; end
        [hfig, hax, sUserData] = ImageCube(ImCube_disp, hdrkwdval, ...
            'fTitleStr', fTitleStr, ...
            'x', plotx, 'y', ploty, 'hfig', hfig);
        
    case 'spread'
        if isempty(hfig), hfig = figure; else, figure(hfig), end
        %[hfig, hax] = PlotSpread;
        [hfig, hax, sUserData] = ImageCubeSpread(ImCube_disp,  hdrkwdval, ...
            'fTitleStr', fTitleStr, ...
            'x', plotx, 'y', ploty, 'hfig', hfig);

    case 'none',
        % no dispaly, do nothing
        
    otherwise,
        error(['unknown plottype ' plottype]);
end

% set all axes to common clim
if isempty(clim),
    clim = AutoClim(ImCube_disp(:), 'one-sided', true);
end
set(hax,'clim',clim)

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

end % main

