function [ImCube, hfig, hax] = FitsPath2ImCube(pn, varargin)
% ImCube = FitsPath2ImCube(pn)
%
% simple routine to collect all the fits files in a subdir
% for example, making an image cube out of all the fits files in a
% phase retrieval subdir

% options
plottype = CheckOption('plottype', 'cube', varargin{:}); % 'spread'
plotx = CheckOption('x', 0, varargin{:});
ploty = CheckOption('y', 0, varargin{:});
xlim = CheckOption('xlim', [], varargin{:});
ylim = CheckOption('ylim', [], varargin{:});

%
listfn = dir(PathTranslator([pn '/*.fits']));
Nf = length(listfn);

% get image size
imtmp = fitsread(PathTranslator([pn '/' listfn(1).name]),'image');
finfo = fitsinfo(PathTranslator([pn '/' listfn(1).name]));
[Nr, Nc] = size(imtmp);

ImCube = zeros([Nf Nr Nc]);
camz = zeros(Nf,1);

ImCube(1,:,:) = imtmp;
camz(1) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'camz');

for ii = 1:Nf

    imtmp = fitsread(PathTranslator([pn '/' listfn(ii).name]),'image');
    finfo = fitsinfo(PathTranslator([pn '/' listfn(ii).name]));

    ImCube(ii,:,:) = imtmp;
    camz(ii) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'camz');
end

switch lower(plottype),
    case 'cube',
        figure,
        [hfig, hax, sUserData] = ImageCube(ImCube, camz,...
            'x', plotx, 'y', ploty);
        
    case 'spread'
        figure
        [hfig, hax] = PlotSpread;
        
    otherwise,
        error(['unknown plottype ' plottype]);
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

